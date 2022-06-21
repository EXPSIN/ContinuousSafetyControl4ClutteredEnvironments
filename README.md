# MPC_counterexapmle

MPC_ContourExample01.m shows that the output of the MPC controller is non-Lipschitz continuous with respect to the system state.

When the mobile agent passes through two obstacles (around t=3.8s), the output of the MPC-based controller changes abruptly.

By modifying the variable Tz in the code, it can be easily found that the non-Lipschitz phenomenon is independent of the sampling interval.

Note: The code has not been tested in versions lower than R2017b.

Tz = 10e-3;
![image](https://user-images.githubusercontent.com/28443522/174713563-62c1eca8-6509-4c0d-8f5d-11ee46b26bb2.png)

Tz = 50e-3;
![image](https://user-images.githubusercontent.com/28443522/174713146-258fc8a3-66f0-4e1a-9977-b9d843fbebf2.png)

Tz = 100e-3;
![image](https://user-images.githubusercontent.com/28443522/174713255-8cd80172-5f0c-45a5-a70c-8dac46cf757e.png)
