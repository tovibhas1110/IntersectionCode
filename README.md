# AiM2XL: Topology Optimization for Wire and Arc Additive Manufacturing (WAAM)

## Overview

This repository contains MATLAB- and Python-based implementations of topology optimization methods developed to generate **additive manufacturing-friendly designs with reduced geometrical complexity**, specifically for the Wire and Arc Additive Manufacturing (WAAM) process.

The work is associated with the **AiM2XL research programme** at Delft University of Technology.

**Research Programme:** [AiM2XL – TU Delft](https://www.tudelft.nl/me/over/afdelingen/materials-science-and-engineering/research/research-programmes/aim2xl)

## Research Motivation

Topology optimization is a powerful computational design method for generating high-performance, lightweight structures. However, conventional topology optimization often produces geometrically complex designs that are difficult to manufacture using WAAM.

WAAM introduces manufacturing considerations related to deposition strategies, geometrical accessibility, and process constraints.

The objective of this research is to develop a topology optimization framework that balances:

- **Structural performance:** Maintaining high stiffness and efficient material distribution.
- **Geometrical simplicity:** Reducing unnecessary geometrical complexity.
- **Manufacturability:** Generating designs better suited to WAAM.
- **Computational efficiency:** Building on established density-based topology optimization methods.

## Methodology

The proposed approach extends the widely used **88-line MATLAB topology optimization code** for compliance minimization.

A key contribution is the use of **stress-field information obtained during compliance minimization to quantify geometrical complexity**.

The methodology involves:

1. Solving the structural compliance minimization problem using density-based topology optimization.
2. Evaluating the stress field associated with the optimized material distribution.
3. Using stress-field information to characterize geometrical complexity.
4. Incorporating complexity-related considerations into the optimization framework.
5. Generating designs that balance structural performance and geometrical simplicity for WAAM applications.

The framework modifies the conventional topology optimization formulation to account for WAAM-related design considerations.

## Implementation

The repository includes MATLAB and Python implementations.

The Python implementation uses numerical and scientific computing libraries for optimization, linear algebra, and visualization.

### Requirements

- Python 3
- NumPy
- SciPy
- Matplotlib

### Getting Started

**1. Clone the repository**

```bash
git clone https://github.com/tovibhas1110/IntersectionCode.git
```

**2. Navigate to the project directory**

```bash
cd IntersectionCode/IntersectionOpt
```

**3. Install the required Python libraries**

```bash
python3 -m pip install numpy scipy matplotlib
```

**4. Run the Python implementation**

```bash
python3 intersection.py
```

## Research Contributions

The primary technical contributions of this work include:

- Extending the classical density-based topology optimization framework to incorporate WAAM-related design considerations.
- Developing a stress-field-based approach for quantifying geometrical complexity.
- Integrating structural performance and geometrical complexity considerations into computational design.
- Implementing numerical optimization algorithms using MATLAB and Python.

## References

**Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011).**

Efficient topology optimization in MATLAB using 88 lines of code.

*Structural and Multidisciplinary Optimization, 43*, 1–16.

https://doi.org/10.1007/s00158-010-0594-7

**AiM2XL Research Programme**

Delft University of Technology, Faculty of Mechanical Engineering.

https://www.tudelft.nl/me/over/afdelingen/materials-science-and-engineering/research/research-programmes/aim2xl

## Research Context

This work contributes to the development of computational design methodologies for large-scale metal additive manufacturing, with particular emphasis on the relationship between structural performance, geometrical complexity, and WAAM manufacturability.

The broader objective is to enable the efficient production of structurally optimized metal components using large-scale additive manufacturing technologies.
