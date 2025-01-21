# M.Sc.-thesis
Code for master's thesis at ETH Zurich, conducted at Harvard Medical School.

Code_minimal contains the minimum working example.

The customized least square algorithm is contained in the Code_minimal/Magnet_solver folder.

least_sqaure_manifold_convergence_cone.m examines how close the initial guess of the magnet's configuration should be to the ground truth to ensure the algorithm's convergence for a given radius of the magnetic sensor circle. The cone is defined on lines 20-21. For example,
'''
theta = linspace(-pi/2+deg2rad(15),pi/2-deg2rad(15),5);
phi = linspace(-pi/2+deg2rad(15),pi/2-deg2rad(15),5);
'''
defines a cone facing upward with a cone size of 75 degrees, the initial guess is the center vector of the cone. Slide 21 shows different sizes of convergence cones with their corresponding sensor configuration. The magnetic sensor configuration is defined on lines 84-89: [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4].

least_sqaure_manifold_simulated_experiment.m examines how the estimation error behaves when there is noise in the measurement. The test configurations are defined on lines 16-71, visualized on slide 25. 210 positions are considered and at each position, 32 orientations are considered. The sensor configuration is defined on lines 84-86. The ground truth magnetic field is computed on lines 123-150, and the artificial noise is added on lines 152-156. Currently, the sensor noise is defined as 5e-8 T on each axis.

Run the program until line 654 shows how much error is contributed from each singular vector. The 4th row in the result, "projected_error/error" shows after projecting the error vector onto each singular vector, how much the projected norm is compared to the error's norm itself. You will always find that the 5th entry is always the largest, which corresponds to the singular vector of the minimum singular value. This is the direction that is the hardest to estimate, which contributes to the most error.

least_sqaure_manifold_physical_experiment.m use the experimental data to estimate the magnet configuration. The experimental data is read on lines 16-116, with the magnet configuration defined, which is the same as on slide 32. The sensor configuration is set on lines 118-139. You can either set a calibrated sensor config (currently it's a calibrated 10cm sensor config) or an uncalibrated one (would be [0.10, 0,0,0,0.10,0,-0.10,0,0,0,-0.10,0]/scale). Run until line 673 will give you the projected error result as in the simulated experiment, which is the result on slide 35. Running until line 723 will give you the mean error for all 32 orientations at each position, which is the result on slide 34. You will observe that even for the same sensor configuration, the estimation error is larger for those configurations whose sigma_min is smaller.
