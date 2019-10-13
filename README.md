# Pressure-equilibrium-phase-composition-for-n-butane-methane

To Run the program simply run the file scripttermproject.m and wait from 6 to 13 seconds to display the results

The objective of the project is to write a Matlab program to predict the behavior of a binary methane and n-butane mixture and construct a plot of pressure vs x1 (methane) at a constant temperature of 100 degrees Fahrenheit. The behavior should be predicted by performing flash calculations of thermodynamic equations obtained from the “Equilibrium constants from a modified Redlich- Kwong equation of state” paper along with the Wilson equation (Wilson, 1968) to predict an initial guess of the equilibrium constants used in the flash calculations.

paper: Soave, G., 1972. Equilibrium constants from a modified Redlich-Kwong equation of state. Chemical engineering science, 27(6), pp.1197-1203.
Link: http://dns2.asia.edu.tw/~ysho/YSHO-English/2000%20Engineering/PDF/Che%20Eng%20Sci27,%201197.pdf

The Wilson equation is the following:
ln(𝐾_i) = ln(𝑃_ci/𝑃) + 5.37(1 + 𝑤_i)(1 − 𝑇_ci/𝑇)

Where,
• 𝐾_i is the ideal equilibrium constant of each component
• 𝑃_ci is the critical pressure of each component
• 𝑇_ci is the critical temperature of each component
• T is the temperature under which the experiment is done
• 𝑊_i is the acentric factor of each component


The following criteria was used to obtain the plot for flash calculations (no modifications were used):
• 𝑍_methane=0.001;𝑍_butane=0.9999 to obtain the first chunk of the graph from the left.
• 𝑍_methane= 0.2; 𝑍_butane=0.80; to obtain the next small chunk of the graph.
• 𝑍_methane= 0.4; 𝑍_butane=0.6 to obtain the near end of the graph.
• 𝐹𝑖𝑛𝑎𝑙𝑦, 𝑍_methane=0.75; 𝑍_butane=0.25 to close the last bit of the graph from the top.
