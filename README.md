# audiocom
**audiocom** is an application written in C for real-time hands-free bidirectional audio communication between two PCs connected over LAN. Audio frames are transmitted using low-latency UDP packets.
<br>It is designed for the development and testing of algorithms for voice communication.

## Features

- **Platform**: Windows
- **Audio Backend**: [PortAudio](http://www.portaudio.com/) (WASAPI interface)
- **FFT Library**: [FFTW](http://www.fftw.org/) for frequency-domain processing
- **Network**: Audio data is sent and received over UDP
- **Configurable**: All parameters are set via a `.txt` configuration file (currently tested only with Sample Rate = 8kHz and Frame size = 128)
- **Audio Enhancement Options**:
  - **AEC (Acoustic Echo Cancellation)**: based on the Multidelay Filter (MDF) [1,2] algorithm
  - **Noise Reduction**: based on log-MMSE [3] spectral estimation.

## Getting Started

1. **Configuration**:
   - Place the `audiocom_095.exe` and the config file `settings_095.txt` in the same directory.
   - Edit the config file to specify:
     - Pre-Shared Key (a unique number to allow authentication at startup),
     - Remote IP address,
     - Local UDP receive port,
     - Remote UDP send port,
   - Ensure the port configuration is complementary on the two PCs.

2. **Network Setup**:
   - Set both network interfaces to **Private Network** mode in Windows.
   - Optionally disable Windows Firewall on both PCs for easier communication.
   - Disable any proprietary audio enhancements (e.g., built-in echo cancellation or noise suppression) in the sound card settings.

3. **Running**:
   - Launch the executable on both PCs after configuring the settings.
   - Type **help** for the list of commands.

## Debugging Tips

If no audio is received or sent:
  - Verify that the PCs are in the same LAN subnet,
  - Check UDP packet transmission using tools like [Wireshark](https://www.wireshark.org/),
  - Confirm correct IP/port settings on both sides,
  - Ensure firewall or antivirus is not blocking the application,
  - Test network connectivity (e.g., via `ping`).

## Future developments

- [x] Acoustic Echo Cancellation (AEC)  
- [x] Noise Cancellation
- [ ] Residual echo suppression 
- [ ] Automatic gain control (AGC)  
- [ ] GUI interface for real-time configuration  
- [ ] Statistics/logging for packet loss and audio quality
- [ ] Integration of forward error correction (FEC) or redundancy schemes to improve reliability
- [ ] Packet Loss Concealment
- [ ] Support for multi-point communication
- [ ] Audio compression to optimize bandwidth consumption.

---

Contributions, bug reports, and feature requests are welcome. Feel free to fork the repository or open an issue!

---
###### [1] - J. S. Soo, K. K. Pang - Multidelay block frequency adaptive filter,  IEEE Transactions on Acoustics Speech Signal Process., Vol. ASSP-38, No. 2, February 1990.
###### [2] - Valin, Jean-Marc - On Adjusting the Learning Rate in Frequency Domain Echo Cancellation With Double-Talk. Audio, Speech, and Language Processing, IEEE Transactions on. 15. 1030 - 1034. 10.1109/TASL.2006.885935 (2007)
###### [3] - Ephraim, Y. and Malah, D. - Speech enhancement using a minimum mean-square error log-spectral amplitude estimator. IEEE Trans. Acoust. Speech, Signal Process., ASSP-23(2), 443-445 (1985)

