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
<img src="/images/HW2_2-DOF_dynamics_estimation_kalman_filter.png" width="50%" height="50%"/>

You can see that every parameters are converge quickly without large overshoot. But also little errors are exist.

## 2.2 Error minimization algorithm
![HW2_2-DOF_dynamics_estimation_minimization_algorithm](/images/HW2_2-DOF_dynamics_estimation_minimization_algorithm.png) | ![HW2_2-DOF_dynamics_estimation_minimization_algorithm2](/images/HW2_2-DOF_dynamics_estimation_minimization_algorithm2.png)
---|---|

Every parameters are converge very quickly than Kalman algorithm and error is almost zero. But in the beginning, very big overshoot exists.

---

# 3.1 1-DOF Joint Space PID CTM Controller
- Target position: 0 deg -> 90 deg
- Target velocity: 30 deg/s
- Use PID controller & Try several gain value

For position controller at critically damped system, P gain K<sub>p</sub> is equal to square of system frequency &omega;<sub>n</sub><sup>2</sup> and D gain K<sub>d</sub> is equal to 2&omega;<sub>n</sub>(&because; &zeta;=1 at critically damped system). If want to decrease steady-state error, we can add I controller.

## 3.1.1 &omega;<sub>n</sub>=5, K<sub>p</sub>=&omega;<sub>n</sub><sup>2</sup>, K<sub>d</sub>=2&omega;<sub>n</sub>, K<sub>i</sub>=0
<img src="/images/HW3_1_wn_5_ki_0.gif" width="50%" height="50%"/><br>
![HW3_1_wn_5_ki_0_pos](/images/HW3_1_wn_5_ki_0_pos.png) | ![HW3_1_wn_5_ki_0_vel](/images/HW3_1_wn_5_ki_0_vel.png)
---|---|

## 3.1.2 &omega;<sub>n</sub>=20, K<sub>p</sub>=&omega;<sub>n</sub><sup>2</sup>, K<sub>d</sub>=2&omega;<sub>n</sub>, K<sub>i</sub>=0
<img src="/images/HW3_1_wn_20_ki_0.gif" width="50%" height="50%"/><br>
![HW3_1_wn_20_ki_0_pos](/images/HW3_1_wn_20_ki_0_pos.png) | ![HW3_1_wn_20_ki_0_vel](/images/HW3_1_wn_20_ki_0_vel.png)
---|---|

## 3.1.3 &omega;<sub>n</sub>=20, K<sub>p</sub>=&omega;<sub>n</sub><sup>2</sup>, K<sub>d</sub>=2&omega;<sub>n</sub>, K<sub>i</sub>=250
<img src="/images/HW3_1_joint_space_1_DOF_CTM_PID_controller.gif" width="50%" height="50%"/><br>
![HW3_1_wn_20_ki_250_pos](/images/HW3_1_joint_space_1_DOF_CTM_PID_controller_pos.png) | ![HW3_1_wn_20_ki_250_vel](/images/HW3_1_joint_space_1_DOF_CTM_PID_controller_vel.png)
---|---|

