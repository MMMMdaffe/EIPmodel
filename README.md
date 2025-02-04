1. Overview
This repository contains MATLAB scripts for simulating the spatiotemporal distributions of ion concentrations and potentials in Electrochemical Ion Pumping (EIP) systems for desalination. These simulations model ion transport and potential distributions in single-electrode and multi-electrode EIP configurations.

2. Installation & Dependencies
Requirements
MATLAB Version: R2023b or later (recommended)
Required Toolboxes: Partial Differential Equation (PDE) Toolbox, Optimization Toolbox
Operating System: Windows/Linux/macOS
No additional dependencies are required beyond standard MATLAB functions.

3. Usage
To run the simulations, download the repository and open MATLAB. Navigate to the folder containing the scripts, then execute the desired .m file using run, for example:
run('Main_single_electrode__CC.m')
Each script generates spatiotemporal data of ion concentrations and potentials.

4. Scripts Description
•	Main_single_electrode__CC.m simulates the spatiotemporal distributions of ion concentrations and electrical potentials in a single-electrode EIP system, where the diluate and concentrate concentrations remain fixed throughout the simulation.

•	Main_single_electrode_desalination.m models the desalination process in a single-electrode EIP system, capturing the dynamic changes in ion concentration as desalination progresses.

•	Main_stack_repeating_circuit.m simulates the spatiotemporal distributions of ion concentrations and electrical potentials in a repeating unit of a multi-electrode EIP stack.

•	Main_stack_terminal.m models the spatiotemporal ion and potential distributions in the terminal circuit of a multi-electrode EIP stack, which differs from repeating units due to electrolysis on the terminal electrode.

6. Please contact the authors if you have any questions. (Weifan Liu, Vanderbilt University, weifan.liu@vanderbilt edu)

