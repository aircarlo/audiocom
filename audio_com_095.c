/*
audio_com_095.c
Version 0.9.5
Changelog: Added a more robust mechanism for connection handling 

Author: Carlo Aironi
Date:   17/04/2025
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <winsock2.h>
#include <windows.h>
#include <portaudio.h>
#include <ws2tcpip.h>
#include <complex.h>
#include <fftw3.h>
#include "a3l_enhancement.h"
#include "a3l_aec.h"

// parameters and default values
int NUM_CHANNELS = 1;
int SAMPLE_RATE = 8000;
int FRAMES_PER_BUFFER = 128;
int IN_PORT = 1000;
int OUT_PORT = 1001;
int BUFFER_SIZE;
int UDP_BUFFER_BYTES = 65536;
int PSK;                            // pre-shared key used for authenticating handshake
int RECV_TIMEOUT = 500;             // timeout in milliseconds to acknowledge that remote host is disconnected
int HANDSHAKE_ATTEMPTS = 15;
char REMOTE_HOST[20] = "0.0.0.0";
int NC_ENABLED = 0;                 // noise cancellation flag
int AEC_ENABLED = 0;                // echo cancellation flag
int AGC_ENABLED = 0;                // automatic gain control flag
int FILTER_LEN;

WSADATA wsaData;
SOCKET inUdpSocket, outUdpSocket;
struct sockaddr_in localAddr, remoteAddr;
PaStream *audioStream;

// struct per il thread audiocom
typedef struct {
    int16_t *sendBuffer;
    int16_t *micBuffer;
    int16_t *aecoutBuffer;
    int16_t *ncoutBuffer;
    int16_t *recvBuffer;
    a3l_aec_state *aec_struct;
    a3l_se_state *nc_struct;
    volatile int *threadRun;
    // FILE *file_1, *file_2, *file_3;                    // DEBUG
} ThreadParams;

// Main Thread I/O + DSP
DWORD WINAPI audiocom_main_thread(LPVOID lpParam) {
    ThreadParams *params = (ThreadParams *)lpParam;
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);																		  

    int16_t *recvBuffer = params->recvBuffer;       // Buffer per i dati ricevuti dal FE
    int16_t *micBuffer = params->micBuffer;         // Buffer per i dati campionati dal MIC
    int16_t *aecoutBuffer = params->aecoutBuffer;   // Buffer per i dati in uscita da AEC
    int16_t *ncoutBuffer = params->ncoutBuffer;     // Buffer per i dati in uscita da NC
    int16_t *sendBuffer = params->sendBuffer;       // Buffer per i dati da inviare al FE
    a3l_aec_state *st_aec = params->aec_struct;
    a3l_se_state *st_nc = params->nc_struct;
    int addrLen = sizeof(remoteAddr);
    int recvCounter = 0;
 
    // Chiamata alle funzioni di inizializzazione delle struct
    a3l_aec_state_init(st_aec, FRAMES_PER_BUFFER, FILTER_LEN, SAMPLE_RATE);
    a3l_se_state_init(st_nc, FRAMES_PER_BUFFER);
    // Dichiarazione dei piani FFTW per AEC
    fftwf_plan fft_plan_aec = fftwf_plan_dft_r2c_1d(st_aec->window_size, st_aec->ifft_buf, st_aec->fft_buf, FFTW_ESTIMATE);
    fftwf_plan ifft_plan_aec = fftwf_plan_dft_c2r_1d(st_aec->window_size, st_aec->fft_buf, st_aec->ifft_buf, FFTW_ESTIMATE);
    // Dichiarazione dei piani FFTW per NC
    fftwf_plan fft_plan_nc = fftwf_plan_dft_r2c_1d(st_nc->window_size, st_nc->ifft_buf, st_nc->fft_buf, FFTW_ESTIMATE);
    fftwf_plan ifft_plan_nc = fftwf_plan_dft_c2r_1d(st_nc->window_size, st_nc->fft_buf, st_nc->ifft_buf, FFTW_ESTIMATE); 
	
    while (*(params->threadRun) && recvCounter<RECV_TIMEOUT) {
        // Lettura di un frame dal microfono
		int framesRecorded = Pa_ReadStream(audioStream, micBuffer, FRAMES_PER_BUFFER);

		// Ricezione di un frame dall'UDP socket di input
        int bytesReceived = recvfrom(inUdpSocket, (char*)recvBuffer, BUFFER_SIZE, 0, NULL, NULL);
        
        // Se è stato ricevuto il frame, lo manda agli altoparlanti
        if (bytesReceived > 0) {
            int framesPlayed = Pa_WriteStream(audioStream, recvBuffer, FRAMES_PER_BUFFER);
            recvCounter = 0;
        } else {
            recvCounter++;
        }

        // Echo Cancellation
        if (AEC_ENABLED) {
            a3l_aec_main_loop(st_aec, recvBuffer, micBuffer, aecoutBuffer, fft_plan_aec, ifft_plan_aec); // AEC
        } else {
            memcpy(aecoutBuffer, micBuffer, FRAMES_PER_BUFFER * sizeof(int16_t)); // PASSTHROUGH
        }

        // Noise Cancellation
        if (NC_ENABLED) {
            a3l_enh_main_loop(st_nc, aecoutBuffer, ncoutBuffer, fft_plan_nc, ifft_plan_nc); // NC
        } else {
            memcpy(ncoutBuffer, aecoutBuffer, FRAMES_PER_BUFFER * sizeof(int16_t)); // PASSTHROUGH
        }

        // AGC (not yet implemented)
        if (AGC_ENABLED) {
            ;
        } else {
            memcpy(sendBuffer, ncoutBuffer, FRAMES_PER_BUFFER * sizeof(int16_t)); // PASSTHROUGH
        }

        // Invio del frame output al socket UDP di uscita
        int bytesSent = sendto(outUdpSocket, (char *)sendBuffer, BUFFER_SIZE, 0, (struct sockaddr *)&remoteAddr, sizeof(remoteAddr));

        // ////////////////// DEBUG ////////////////////
        // // scrittura su file(s) bin
        // float frame2file_1[128], frame2file_2[128], frame2file_3[128];
        // for (int i = 0; i < 128; i++) {
        //     frame2file_1[i] = micBuffer[i] / 32768.0f;
        //     frame2file_2[i] = recvBuffer[i] / 32768.0f;
        //     frame2file_3[i] = sendBuffer[i] / 32768.0f;
        // }
        // fwrite(frame2file_1, sizeof(float), 128, params->file_1);
        // fwrite(frame2file_2, sizeof(float), 128, params->file_2);
        // fwrite(frame2file_3, sizeof(float), 128, params->file_3);
        // /////////////////////////////////////////////

    }
    if (recvCounter==RECV_TIMEOUT) {
        printf("Communication closed by remote host.\n");
        printf(">> ");
    }
    return 0;
}

// Lettura dei parametri dal file di configurazione
void parseParams(const char *filename, int *psk, char* remoteHost, int *inPort, int *outPort, int *sampleRate, int *framesPerBuffer, int *UDP_buffer_bytes, int *nc_enabled, int *aec_enabled, int *aec_filter_len){
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Errore nell'apertura del file %s\n", filename);
        exit(1);
    }

    char line[256];
    while (fgets(line, sizeof(line), file))
    {
        line[strcspn(line, "\n")] = '\0';
        
        if (strncmp(line, "PSK", 3) == 0){
            sscanf(line, "PSK %d", psk);
        }
        else if (strncmp(line, "REMOTE_HOST", 11) == 0){
            sscanf(line, "REMOTE_HOST %s", remoteHost);
        }
        else if (strncmp(line, "IN_PORT", 7) == 0){
            sscanf(line, "IN_PORT %d", inPort);
        }
        else if (strncmp(line, "OUT_PORT", 8) == 0){
            sscanf(line, "OUT_PORT %d", outPort);
        }
        else if (strncmp(line, "SAMPLE_RATE", 11) == 0){
            sscanf(line, "SAMPLE_RATE %d", sampleRate);
        }
        else if (strncmp(line, "FRAMES_PER_BUFFER", 17) == 0){
            sscanf(line, "FRAMES_PER_BUFFER %d", framesPerBuffer);
        }
        else if (strncmp(line, "UDP_BUFFER_BYTES", 16) == 0){
            sscanf(line, "UDP_BUFFER_BYTES %d", UDP_buffer_bytes);
        }
        else if (strncmp(line, "NC", 2) == 0){
            sscanf(line, "NC %d", nc_enabled);
        }
        else if (strncmp(line, "AEC", 3) == 0){
            sscanf(line, "AEC %d", aec_enabled);
        }
        else if (strncmp(line, "FILTER_LEN", 10) == 0){
            sscanf(line, "FILTER_LEN %d", aec_filter_len);
        }
    }
    fclose(file);
}

// Funzione per la lettura dell'indirizzo IP delle interfacce NIC
void printIP(){
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    struct addrinfo hints = {0}, *res, *p;
    hints.ai_family = AF_INET; // Solo IPv4
    getaddrinfo(hostname, NULL, &hints, &res);

    printf("- Local host IP:\t\t");
    for (p = res; p != NULL; p = p->ai_next) {
        struct sockaddr_in *addr = (struct sockaddr_in *)p->ai_addr;
        printf("%s   ", inet_ntoa(addr->sin_addr));
    }
    printf("\n");

    freeaddrinfo(res);
}

// Print startup info
void printStartup() {
    printf("      _                                                             \n");
    printf("  /\\  _) |   _. |_     _.      _| o  _   _  _  ._ _     _. ._  ._  \n");
    printf(" /--\\ _) |_ (_| |_)   (_| |_| (_| | (_) (_ (_) | | |   (_| |_) |_) \n");
    printf("                                                           |   |    \n");
    printf("                                                                    \n");
    printf("====================================================================\n");
    printf("Version 0.9.5                                   Author: Aironi Carlo\n");
    printf("====================================================================\n");
    printf("\n");
    printIP();
    printf("- Remote host IP:\t\t%s\n", REMOTE_HOST);
    printf("- Input UDP port:\t\t%d\n", IN_PORT);
    printf("- Output UDP port:\t\t%d\n", OUT_PORT);
    printf("- Audio sample rate:\t\t%d Hz\n", SAMPLE_RATE);
    printf("- Audio frame size:\t\t%d\n", FRAMES_PER_BUFFER);
    printf("- UDP input socket buffer size: %d bytes\n", UDP_BUFFER_BYTES);
    if (NC_ENABLED){
        printf("- Noise Cancellation:\t\tON\n");
    } else {
        printf("- Noise Cancellation:\t\tOFF\n");
    }
    if (AEC_ENABLED){
        printf("- Echo Cancellation:\t\tON - ");
        printf("AEC filter length: %d\n", FILTER_LEN);
    } else {
        printf("- Echo Cancellation:\t\tOFF\n");
    }
    printf("\n");
}

void printHelp() {
    printf("\n");
    printf("The program reads the startup parameters from settings_095.txt, located in the same directory as audio_com_095.exe.\n");
    printf("Make sure you have entered the arguments correctly before starting the program:\n");
    printf("- PSK\t\t\t(int)\tPre-Shared Key for authentication, must be the same on both terminals\n");
    printf("- REMOTE_HOST\t\t(str)\tIP address of the remote host in the format [XXX.XXX.XXX.XXX]\n");
    printf("- IN_PORT\t\t(int)\tInput UDP port, must match the OUT_PORT of remote host\n");
    printf("- OUT_PORT\t\t(int)\tOutput UDP port, must match the IN_PORT of remote host\n");
    printf("- SAMPLE_RATE\t\t(int)\tAudio sample rate (default 8000)\n");
    printf("- FRAMES_PER_BUFFER\t(int)\tAudio frame size in samples (default 128)\n");
    printf("- UDP_BUFFER_BYTES\t(int)\tWinsock internal UDP buffer (default 1024, larger values may introduce latency)\n");
    printf("- NC\t\t\t(int)\tFlag for activating [1] or deactivating [0] Noise Cancellation\n");
    printf("- AEC\t\t\t(int)\tFlag for activating [1] or deactivating [0] Echo Cancellation\n");
    printf("- FILTER_LEN\t\t(int)\tAEC adaptive filter length (must be a power of 2)\n");
    printf("\n");
    printf("Commands:\n");
    printf("- 'start'\t\tStarts communication\n");
    printf("- 'stop\t\t\tEnds communication\n");
    printf("- 'exit'\t\tClose the program\n");
    printf("- 'help'\t\tShows help\n");
    printf("\n");
}

int main() {

    // Read parameters from file
    parseParams("settings_095.txt",
                &PSK, 
                REMOTE_HOST,
                &IN_PORT,
                &OUT_PORT,
                &SAMPLE_RATE,
                &FRAMES_PER_BUFFER,
                &UDP_BUFFER_BYTES,
                &NC_ENABLED,
                &AEC_ENABLED,
                &FILTER_LEN);

    BUFFER_SIZE = (FRAMES_PER_BUFFER * sizeof(short));

    WSAStartup(MAKEWORD(2, 2), &wsaData);
   
    printStartup();

    u_long mode = 1; // modalità dei socket UDP (1 = non-blocking)

    // Creazione socket UDP di ricezione
    inUdpSocket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    ioctlsocket(inUdpSocket, FIONBIO, &mode);
    localAddr.sin_family = AF_INET;
    localAddr.sin_port = htons(IN_PORT);
    localAddr.sin_addr.s_addr = INADDR_ANY;
    bind(inUdpSocket, (struct sockaddr*)&localAddr, sizeof(localAddr));
    setsockopt(inUdpSocket, SOL_SOCKET, SO_RCVBUF, (char*)&UDP_BUFFER_BYTES, sizeof(UDP_BUFFER_BYTES));

    // Creazione socket UDP di trasmissione
    outUdpSocket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    ioctlsocket(outUdpSocket, FIONBIO, &mode);
    remoteAddr.sin_family = AF_INET;
    remoteAddr.sin_port = htons(OUT_PORT);
    remoteAddr.sin_addr.s_addr = inet_addr(REMOTE_HOST);
    setsockopt(outUdpSocket, SOL_SOCKET, SO_SNDBUF, (char*)&UDP_BUFFER_BYTES, sizeof(UDP_BUFFER_BYTES));

    // Allocazione e inizializzazione dei buffer audio
    int16_t *sendBuffer = (int16_t *)malloc(FRAMES_PER_BUFFER * sizeof(int16_t));
    int16_t *micBuffer = (int16_t *)malloc(FRAMES_PER_BUFFER * sizeof(int16_t));
    int16_t *aecoutBuffer = (int16_t *)malloc(FRAMES_PER_BUFFER * sizeof(int16_t));
    int16_t *ncoutBuffer = (int16_t *)malloc(FRAMES_PER_BUFFER * sizeof(int16_t));
    int16_t *recvBuffer = (int16_t *)malloc(FRAMES_PER_BUFFER * sizeof(int16_t));
    memset(sendBuffer, 0, FRAMES_PER_BUFFER * sizeof(int16_t));
    memset(micBuffer, 0, FRAMES_PER_BUFFER * sizeof(int16_t));
    memset(aecoutBuffer, 0, FRAMES_PER_BUFFER * sizeof(int16_t));
    memset(ncoutBuffer, 0, FRAMES_PER_BUFFER * sizeof(int16_t));
    memset(recvBuffer, 0, FRAMES_PER_BUFFER * sizeof(int16_t));

    // ////////////////// DEBUG ////////////////////
    // // Apertura del file per l'output su .bin
    // FILE *file_1 = fopen("signal_s_RAW.bin", "ab");
    // FILE *file_2 = fopen("signal_x_RAW.bin", "ab");
    // FILE *file_3 = fopen("signal_e_RAW.bin", "ab");
    // /////////////////////////////////////////////

    a3l_aec_state *st_aec = (a3l_aec_state*)malloc(sizeof(a3l_aec_state));
    a3l_se_state *st_nc = (a3l_se_state*)malloc(sizeof(a3l_se_state));

    // Inizializzazione della struttura per i parametri del thread audiocom
    volatile int threadRun;     // variabile di controllo per il ciclo while nel thread
    ThreadParams *io_params = malloc(sizeof(ThreadParams));
    io_params->sendBuffer = sendBuffer;
    io_params->micBuffer = micBuffer;
    io_params->aecoutBuffer = aecoutBuffer;
    io_params->ncoutBuffer = ncoutBuffer;
    io_params->recvBuffer = recvBuffer;
    io_params->aec_struct = st_aec;
    io_params->nc_struct = st_nc;
    io_params->threadRun = &threadRun;
    // io_params->file_1 = file_1;                      // DEBUG
    // io_params->file_2 = file_2;                      // DEBUG
    // io_params->file_3 = file_3;                      // DEBUG

    // wait for user prompt
    char userInput[256];
    printf(">> Type 'help' for the list of commands\n");
prompt:     // goto flag, NON ABUSARNE!
    printf(">> ");
    // Legge l'input dell'utente (fino a 255 caratteri + '\0')
    fgets(userInput, sizeof(userInput), stdin);
    // Rimuove il carattere di newline alla fine della stringa
    userInput[strcspn(userInput, "\n")] = '\0';
    // Verifica se l'input è 'start'
    if (strcmp(userInput, "start") == 0) {
        // esce dall'IF e avvia il thread di comunicazione
        threadRun = 1;
    } else if (strcmp(userInput, "stop") == 0) {
        printf(">> Communication closed\n");
        threadRun = 0;
        goto prompt;
    } else if (strcmp(userInput, "exit") == 0) {
        printf(">> BYE\n");
        Sleep(500);
        goto theend;
    } else if (strcmp(userInput, "help") == 0) {
        printHelp();
        goto prompt;
    } else {
        printf(">> Invalid command\n");
        goto prompt;
    }

    // Start connection...
    printf(">> Awaiting response from remote host...");
    // clear Winsock UDP cache (possibly useless)
    for (int i=0; i<10; i++){
        recvfrom(inUdpSocket, (char*)recvBuffer, BUFFER_SIZE, 0, NULL, NULL);
    }
    // Send handshake request (a frame with the Pre-Shared Key) and wait for handshake response
    int handshakeReceived = -1;
    int handshakeAttempts = 0;
    sendBuffer[0] = PSK;
    while(handshakeReceived == -1 && handshakeAttempts < HANDSHAKE_ATTEMPTS) {
        Sleep(1000);
        int bytesSent = sendto(outUdpSocket, (char *)sendBuffer, BUFFER_SIZE, 0, (struct sockaddr *)&remoteAddr, sizeof(remoteAddr));
        handshakeReceived = recvfrom(inUdpSocket, (char*)recvBuffer, BUFFER_SIZE, 0, NULL, NULL);
        handshakeAttempts++;
    }

    // exit after HANDSHAKE_ATTEMPTS unsuccessful attempts
    if (handshakeReceived == -1) {
        printf(" fail (timeout)\n");
        goto prompt;
    }

    if (recvBuffer[0] != PSK) {
        printf(" fail (keys mismatch)\n");
        goto prompt;
    } else {
        printf(" ok\n");
        printf(">> Communication started, ");
    }

    // Inizializzazione PortAudio
    Pa_Initialize();
    Pa_OpenDefaultStream(&audioStream, 1, 1, paInt16, SAMPLE_RATE, FRAMES_PER_BUFFER, NULL, NULL);
    Pa_StartStream(audioStream);
    // Avvia il thread
    HANDLE comThread = CreateThread(NULL, 0, audiocom_main_thread, io_params, 0, NULL);
    printf("type 'stop' to terminate\n");
    goto prompt;
    WaitForSingleObject(comThread, INFINITE);

theend:    
    // Cleanup
    Pa_StopStream(audioStream);
    Pa_CloseStream(audioStream);
    Pa_Terminate();
    closesocket(inUdpSocket);
    closesocket(outUdpSocket);
    WSACleanup();
    a3l_aec_state_free(st_aec);
    a3l_se_state_free(st_nc);
    free(sendBuffer);
    free(micBuffer);
    free(aecoutBuffer);
    free(ncoutBuffer);
    free(recvBuffer);

    system("pause");    // evita la chiusura automatica della finestra
    return 0;
}