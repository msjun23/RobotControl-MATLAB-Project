# RobotControl-MATLAB-Project

**Kwangwoon University, Seoul, South Korea**
<br>
Term project for Robot Control(Second semester, 2021.)

PDF, lecture material is own Kwangwoon Univ.

---

# 1. 3-DOF links free fall simulation
Code is [here](https://github.com/msjun23/RobotControl-MATLAB-Project/blob/main/RobotControl_TermProject/HW1_3_DOF_simulation.m) and dynamics modeling using Lagrangian function is [here](https://github.com/msjun23/RobotControl-MATLAB-Project/blob/main/RobotControl_TermProject/HW1_3_DOF_dynamics.m).

<img src="/images/HW1_3_DOF_simulation.gif" width="50%" height="50%"/>

---

# 2. 2-DOF Dynamics parameter estimation
Code is [here](https://github.com/msjun23/RobotControl-MATLAB-Project/blob/main/RobotControl_TermProject/HW2_2_DOF_dynamics_parameter_estimation.m) and two link dynamics equation is [here](https://github.com/msjun23/RobotControl-MATLAB-Project/blob/main/RobotControl_TermProject/HW2_two_link.m).

To estimation, two algorithms are available.

## 2.1 Kalman Filter based parameter estimation algorithm
<img src="/images/2-DOF_dynamics_estimation_kalman_filter.png" width="50%" height="50%"/>

You can see that every parameters are converge quickly without large overshoot. But also little errors are exist.

## 2.2 Error minimization algorithm
![2-DOF_dynamics_estimation_minimization_algorithm](/images/2-DOF_dynamics_estimation_minimization_algorithm.png) | ![2-DOF_dynamics_estimation_minimization_algorithm2](/images/2-DOF_dynamics_estimation_minimization_algorithm2.png)
---|---|

Every parameters are converge very quickly than Kalman algorithm and error is almost zero. But in the beginning, very big overshoot exists.
