# Rice-Fading-Matlab
The Rician fading model implemented in Matlab and was created using Matlab 2018a.
Plots the theoretical and simulated, envelope and phase porbability density functions (PDFs)

The Rician fading model describes the small-scale fading ocurring between and transmitter and receiver antenna pair, where a scattered and dominant signal components exist.
A special case of this model is Rayleigh fading, when only a scattered signal component exists. This is represented by K = 0.

Run main.m to start the GUI if Matlab is already installed.
Alternatively if Matlab isn't installed, can run the installer from the build folder, which requires an internet connection to download the required files.

When running the program the intial window appears:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Rician-Fading-Matlab/main/docs/window.png)

Entering values for the Rician K factor, the root mean sqaure of the signal (\hat{r}^2), and \phi the phase parameter:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Rician-Fading-Matlab/main/docs/inputs.png)

The theoretical evenlope PDF and phase PDFs are plotted to compare with the simulation and gives the execution time for the theoretical calculations and simulations together:

![ScreenShot](https://raw.github.com/Jonathan-Browning/Rician-Fading-Matlab/main/docs/results.png)
