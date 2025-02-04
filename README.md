# Pressure-equilibrium-phase-composition-for-n-butane-methane

To run the program, simply execute the scripttermproject.m file and wait 6 to 13 seconds for the results to display.

This MATLAB program is designed to predict the behavior of a binary mixture of methane and n-butane. The program constructs a plot of pressure versus the mole fraction of methane (x1) at a constant temperature of 100°F. The behavior is predicted using flash calculations derived from thermodynamic equations based on the paper "Equilibrium Constants from a Modified Redlich-Kwong Equation of State" by Soave (1972), alongside the Wilson equation (Wilson, 1968), to estimate the equilibrium constants for the flash calculations.

Reference Paper: Soave, G., 1972. Equilibrium constants from a modified Redlich-Kwong equation of state. Chemical engineering science, 27(6), pp.1197-1203.
Link: http://dns2.asia.edu.tw/~ysho/YSHO-English/2000%20Engineering/PDF/Che%20Eng%20Sci27,%201197.pdf

Wilson Equation:
ln(𝐾_i) = ln(𝑃_ci/𝑃) + 5.37(1 + 𝑤_i)(1 − 𝑇_ci/𝑇)

Where,
• 𝐾_i is the ideal equilibrium constant of each component
• 𝑃_ci is the critical pressure of each component
• 𝑇_ci is the critical temperature of each component
• T is the temperature under which the experiment is done
• 𝑊_i is the acentric factor of each component


Flash Calculation Criteria:
The following values were used to generate the plot:
• 𝑍_methane=0.001;𝑍_butane=0.9999 to obtain the first chunk of the graph from the left.
• 𝑍_methane= 0.2; 𝑍_butane=0.80; to obtain the next small chunk of the graph.
• 𝑍_methane= 0.4; 𝑍_butane=0.6 to obtain the near end of the graph.
• 𝐹𝑖𝑛𝑎𝑙𝑦, 𝑍_methane=0.75; 𝑍_butane=0.25 to close the last bit of the graph from the top.
