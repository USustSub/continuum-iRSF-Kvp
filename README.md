# continuum-iRSF-Kvp

Continuum plasticity model for Sequences of Earthquakes and Aseismic Slip (SEAS) using an invariant, consistently-linearized formulation of rate-and-state friction. This models also features Kelvin-viscoplasticity for regularization of the fault width. However, obtaining seismic slip rates then requires grid sizes down to meters. When the Kelvin viscosity is set to 0 that regularisation is no active. However, strain-rate dependent properties are already grid-size independent, because of the rate-dependence introduced by the calculation of the slip rate ( = plastic strain rate * fault width of 2 * grid size ). This model is compared to a discrete Finite Element Method model, which reproduces is results successfully. Explanations and comparison are available in Gourdarzi et al., Computational Mechanics, 2023. 

Goudarzi, M., Gerya, T. & Dinther, Y. van. A comparative analysis of continuum plasticity, viscoplasticity and phase-field models for earthquake sequence modeling. Comput. Mech. 1–19 (2023) doi:10.1007/s00466-023-02311-0.
See https://link.springer.com/article/10.1007/s00466-023-02311-0

This code has been developed as part of the research programme DeepNL, financed by the Dutch Research Council (NWO) under project number DeepNL.2018.033.


