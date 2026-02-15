# **LLG_Solver_CPP — Landau–Lifshitz–Gilbert Equation Solver**

An authentic, modular C++ framework for solving the Landau–Lifshitz–Gilbert (LLG) equation. This project is designed for micromagnetic simulations, focusing on numerical stability, clean architecture, and extensibility.

---

## **1. Introduction**

The dynamics of magnetization $\mathbf{M}$ are governed by the **Landau–Lifshitz–Gilbert (LLG)** equation:

$$\frac{d\mathbf{M}}{dt} = -\gamma (\mathbf{M} \times \mathbf{H}_{\text{eff}}) + \frac{\alpha}{M_s} \left( \mathbf{M} \times \frac{d\mathbf{M}}{dt} \right)$$

Where:
* $\gamma$ — Gyromagnetic ratio
* $\alpha$ — Gilbert damping parameter
* $M_s$ — Saturation magnetization
* $\mathbf{H}_{\text{eff}}$ — Effective magnetic field (External, Exchange, Anisotropy, etc.)



---

## **2. Features**

* **Modern C++17:** Structured with clear header/source separation in `include/` and `src/`.
* **Modular Fields:** Independent modules for Zeeman, Exchange, and Anisotropy fields.
* **Numerical Integrators:** Implementation of **RK4** (Runge-Kutta 4th Order) for high precision.
* **Normalization Enforcement:** Automatically maintains $|\mathbf{M}| = M_s$ after each time step to ensure physical validity.
* **Hysteresis Modeling:** Built-in capability to simulate magnetization loops at varying temperatures.

---

## **3. Repository Structure**



```text

LLG_Solver_CPP/
├── analysis
│   └── plot_magnetization.py
├── apps
│   └── run_hysteresis.cpp
├── build
│   └── a.out
├── data
├── include
│   ├── DelmLlg.h
│   ├── ExchangeField.h
│   ├── ExternalField.h
│   ├── LlgSolver.h
│   └── utils.h
├── README.md
├── run.sh
├── scripts
├── src
│   ├── DelmLlg.cpp
│   ├── ExchangeField.cpp
│   ├── LlgSolver.cpp
│   ├── main.cpp
│   └── utils.cpp
└── tests
    ├── hystersis_vr2.cpp
    ├── LLG_single_spin_T_0Kelvin.cpp
    ├── LLG_single_spin_T_finite.cpp
    └── LLG_solver_2D_grid.cpp
