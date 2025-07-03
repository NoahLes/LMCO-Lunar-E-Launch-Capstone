Arduino Code for Subscale demonstration of E-Launcher:

LMCO_FIRE_CONTROL - Utilizes infared LED's to track projectile across the launcher barrel. 
Calculates time and velocity data of the projectile in order to dynamically switch on each stage 
of the launcher to increase the forward force from the magnetic field produced within each stage.

LMCO_DAQ_CONTROLPANEL - Produces a web interface which allows the user to wirelessly display Launcher 
diagnostics via laptop or mobile phone. Users are able to download diagnostics in .csv format. "FIRE" 
button implemented which allows the user to fire the launcher without physically interacting with the launcher.

MATLAB Simulation Code for the Conceptual and Subscale variants of our E-Launcher:

multistage_conceptual - MATLAB script that simulates the magnetic field effects and other physical aspects of
the full scale lunar integration of the launcher. Does not take into account factors such as thermal properties, 
and the vaccum of space. Provides us with simulated diagnostics of the launcher such as capacitor voltage curves, 
exit velocity of the projectile, current through the coils, and magnetic field properties.

Supporting Files: ConceptualBzField, ConceptualBrField

multistage_subscale - MATLAB script that simulates the magnetic field effects and other physical aspects of the subscale
integration of the launcher. Does not take intop account factors such as thermal properties, and the air drag on the
projectile. Provides us with simulated diagnostics of the launcher such as capactior voltage curves, exit velocity of the
projectile, current through the coils, and magnetic field properties. Simulation can be validated by comparing simulation output
values with the output diagnostics of the subscale demonstration.

Supporting Files: SubscaleBzField, SubscaleBrField, INPUTPARAMETERS
