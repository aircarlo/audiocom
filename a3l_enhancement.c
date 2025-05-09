/*
File:       a3l_enhancement.c
Version:    2.0
Author:     Aironi Carlo
Date:       07feb 2025

Implementation of the logMMSE algorithm [1].
Uses the exponential integral piecewise approximation (Martin et al. [2][3])
Takes into account Speech Presence Uncertainty

[1] Ephraim, Y. and Malah, D. (1985). Speech enhancement using a minimum 
    mean-square error log-spectral amplitude estimator. IEEE Trans. Acoust., 
    Speech, Signal Process., ASSP-23(2), 443-445.
[2] R. Martin and D. Malah and R. V. Cox and A. J. Accardi, "A Noise 
    Reduction Preprocessor for Mobile Voice Communication,” Technion - Israel
    Institute of Technology, Technical Report, CCIT No. 459, Dec. 2003.
[3] Ephraim, Y. and Cohen, I. (2006), Recent advancements in speech 
    enhancement, in Dorf, R. C. (Ed.), The Electrical Engineering Handbook, 
    Boca Raton, FL: CRC Press.

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "a3l_enhancement.h"

#ifndef M_PI
#define M_PI 3.141592f
#endif

// Parametri del cancellatore di rumore
#define ALPHA 0.98f                             // Coefficiente per SNR a priori
#define CSI_MIN 0.0032f          				// Limite minimo per csi_k (-25dB), preso dall'implementazione di Loizou, ma non in Ephraim-Malah
#define ETA 0.15f               				// Soglia VAD
#define MU 0.98f               					// Fattore di aggiornamento del noise_PSD
#define SPU 0                                   // Speech Presence Uncertainty flag
#define QKR 2.3333f                             // SPU parameter


// Allocazione dinamica e inizializzazione della struttura
void a3l_se_state_init(a3l_se_state *st, int frame_size) {
	st->frame_size = frame_size;
	st->window_size = 2 * frame_size;
	st->vad_decision = 0.0f;
    st->noise_PSD_buf = (float*)calloc(st->window_size, sizeof(float));
	for (int i = 0; i < st->window_size; i++) {
		st->noise_PSD_buf[i] = 0.1f;
	}
	st->Y_psd_buf = (float*)calloc(st->window_size, sizeof(float));
	st->yf = (float*)calloc(st->window_size, sizeof(float));
    st->yf_buf = (float*)calloc(st->frame_size, sizeof(float));
    st->out_buf = (float*)calloc(st->frame_size, sizeof(float));
	st->xf_hat = (float*)calloc(st->window_size, sizeof(float));
    st->hann_win = (float*)calloc(st->window_size, sizeof(float));
	st->Yc = (float _Complex*)calloc(st->window_size, sizeof(float _Complex));
	st->Y_mag = (float*)calloc(st->window_size, sizeof(float));
	st->Y_psd = (float*)calloc(st->window_size, sizeof(float));
	st->gamma_k = (float*)calloc(st->window_size, sizeof(float));
	st->csi_k = (float*)calloc(st->window_size, sizeof(float));
	st->log_sigma_k = (float*)calloc(st->window_size, sizeof(float));
	st->Y_pha = (float*)calloc(st->window_size, sizeof(float));
	st->X_hat = (float _Complex*)calloc(st->window_size, sizeof(float _Complex));
	st->A = (float*)calloc(st->window_size, sizeof(float));
	st->vk = (float*)calloc(st->window_size, sizeof(float));
	st->ei_vk = (float*)calloc(st->window_size, sizeof(float));
	st->hw = (float*)calloc(st->window_size, sizeof(float));
    st->ifft_buf = (float*)calloc(st->window_size, sizeof(float));
    st->fft_buf = (float _Complex*)calloc(st->window_size, sizeof(float _Complex));
    // inizializzazione della finestra con il profilo di Hanning
    for (int i = 0; i < st->window_size; i++) {
        st->hann_win[i] = 0.5f - 0.5f * cosf(2.0f * M_PI * ((float)i / (float)(st->window_size - 1)));
    }
}

// Deallocazione della struttura
void a3l_se_state_free(a3l_se_state *st) {
    // Dealloca la memoria per ciascun array membro
	free(st->noise_PSD_buf);
	free(st->Y_psd_buf);
	free(st->yf);
    free(st->yf_buf);
    free(st->out_buf);
	free(st->xf_hat);
	free(st->Yc);
	free(st->Y_mag);
	free(st->Y_psd);
	free(st->gamma_k);
	free(st->csi_k);
	free(st->log_sigma_k);
	free(st->Y_pha);
	free(st->X_hat);
    free(st->A);
    free(st->vk);
    free(st->ei_vk);
    free(st->hw);
    free(st->fft_buf);
    free(st->ifft_buf);
    free(st);
}

// Funzione che rende coniugato simmetrico un array complesso di lunghezza L (pari)
static void conj_sym(float _Complex *Arr, int L) {
    // L è la lunghezza totale dell'array (uguale alla lunghezza del segnale originale)
    // Le prime 129 posizioni contengono i campioni dello spettro (frequenze positive)
    // Copia e coniuga la parte del secondo spettro (frequenze negative)
    for (int i = 1; i <= L/2; i++) {
        Arr[L-i] = conjf(Arr[i]);
    }
}

// Funzione che moltiplica per una costante tutti gli elementi di un array a valori reali di lunghezza L
// K1 = numeratore
// K2 = denominatore
static void arr_scaling_re(float *Arr, int L, float K1, float K2) {
    for (int i = 0; i < L; i++) {
        Arr[i] *= K1;
        Arr[i] /= K2;
    }
}

// Loop di cancellazione del rumore
void a3l_enh_main_loop(a3l_se_state *st, const int16_t *in, int16_t *out, fftwf_plan fft_plan, fftwf_plan ifft_plan) {

    for (int i = 0; i < st->frame_size; i++) {
        st->yf[i] = st->yf_buf[i];
        st->yf[i + st->frame_size] = (float)in[i] / 32768.0f;
        st->yf_buf[i] = (float)in[i] / 32768.0f;                   // salva nel buffer per il ciclo successivo
    }

    for (int i = 0; i < st->window_size; i++) {
        st->yf[i] = st->yf[i] * st->hann_win[i];
    } 

    // FFT
    memcpy(st->ifft_buf, st->yf, st->window_size * sizeof(float));
    fftwf_execute(fft_plan);                                                        // FFT
    memcpy(st->Yc, st->fft_buf, st->window_size * sizeof(float _Complex));
    conj_sym(st->Yc, st->window_size);                                            // rendo l'output coniugato simmetrico, poichè il piano fftw non lo fa
    
    st->vad_decision = 0.0f;
    for (int i = 0; i < st->window_size; i++) {
        // Split real and imaginary parts
        float Y_re = crealf(st->Yc[i]);
        float Y_im = cimagf(st->Yc[i]);
        // Compute magnitude and phase
        st->Y_mag[i] = cabsf(st->Yc[i]);
        st->Y_psd[i] = powf(st->Y_mag[i], 2.0f);
        st->Y_pha[i] = atan2f(Y_im, Y_re);
        st->gamma_k[i] = fminf(st->Y_psd[i] / st->noise_PSD_buf[i], 40.0f);		// gamma_k
        // D-D approach to estimate a priori SNR
        st->csi_k[i] = ALPHA * st->Y_psd_buf[i] / st->noise_PSD_buf[i] + (1 - ALPHA) * fmaxf(st->gamma_k[i] - 1.0f, 0.0f);
        st->csi_k[i] = fmaxf(CSI_MIN, st->csi_k[i]);
        // VAD
        st->log_sigma_k[i] = st->gamma_k[i] * st->csi_k[i] / (1.0f + st->csi_k[i]) - logf(fmaxf(1.0f + st->csi_k[i], 1e-6f));
        st->vad_decision += (st->log_sigma_k[i] / st->window_size);
    }

    // Update the PSD noise estimate if vad_decision is below the threshold
    if (st->vad_decision < ETA) {
        for (int i = 0; i < st->window_size; i++) {
            st->noise_PSD_buf[i] = MU * st->noise_PSD_buf[i] + (1.0f - MU) * st->Y_psd[i];
        }
    }

    for (int i = 0; i < st->window_size; i++) {
        // Compute gain function
        st->A[i] = st->csi_k[i] / (1.0f + st->csi_k[i]);
        st->vk[i] = st->A[i] * st->gamma_k[i];
        // Exponential integral calculation function with piecewise approximation (Martin et al.)
        if (st->vk[i] < 0.1f) {
            st->ei_vk[i] = -2.31f * log10f(st->vk[i]) - 0.6f;
        } else if (st->vk[i] >= 0.1f && st->vk[i] <= 1.0f) {
            st->ei_vk[i] = -1.544f * log10f(st->vk[i]) + 0.166f;
        } else {
            st->ei_vk[i] = powf(10.0f, (-0.52f * st->vk[i]) - 0.26f);
        }

        st->hw[i] = st->A[i] * expf((0.5f * st->ei_vk[i]));
        // Gain application
#if SPU
        float L = QKR * expf(st->vk[i]) / (1.0f + st->csi_k[i]);
        st->Y_mag[i] = st->Y_mag[i] * st->hw[i] * L / (1.0f + L);
#else
        st->Y_mag[i] *= st->hw[i];
#endif
        // Save for a-priori SNR estimation in next frame
        st->Y_psd_buf[i] = (st->Y_mag[i] * st->Y_mag[i]);
        // Complex spectrum reassembly with filtered Y_mag and noisy Y_pha
        st->X_hat[i] = st->Y_mag[i] * cexpf(I * st->Y_pha[i]);          // Evita calcoli separati di seno e coseno

    }

    // Inverse FFT
    memcpy(st->fft_buf, st->X_hat, st->window_size * sizeof(float _Complex));
    fftwf_execute(ifft_plan);
    memcpy(st->xf_hat, st->ifft_buf, st->window_size * sizeof(float));
    arr_scaling_re(st->xf_hat, st->window_size, 1.0f, (float)st->window_size);

    // re-assembling output
    for (int i = 0; i < st->window_size; i++) {											// clipping and scaling (float)
        if (st->xf_hat[i] > 1.0f){
            st->xf_hat[i] = 1.0f;
        } else if(st->xf_hat[i] < -1.0f){
            st->xf_hat[i] = -1.0f;
        }
    }

    for (int i = 0; i < st->frame_size; i++) {											// cast to short
        out[i] = (int16_t)((st->xf_hat[i] + st->out_buf[i]) * 32768.0f);    			// OLA
        st->out_buf[i] = st->xf_hat[i + st->frame_size];								// save output into buffer
    }

}
