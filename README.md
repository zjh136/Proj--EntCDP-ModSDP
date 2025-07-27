# EntCDP-ModSDP includes two complementary models

- **EntCDP** is designed to identify **common signaling pathways** across multiple cancer types.  
  ➤ To run EntCDP, execute:  
  `run_EntCDP.m`

- **ModSDP** aims to identify **specific signaling pathways** in one group of cancers **relative to another group**.  
  ➤ To run ModSDP, execute:  
  `run_ModSDP.m`

---

## Requirements

- **MATLAB** (tested on R2021a and above)
- **IBM ILOG CPLEX Optimizer** (for solving ModSDP via mathematical programming)
- Alternatively, we provide a Genetic Algorithm (GA)-based solver (ModSDP_GA_matlab.m) in this package as a substitute for CPLEX.
  [Download CPLEX from IBM](https://www.ibm.com/products/ilog-cplex-optimization-studio)
