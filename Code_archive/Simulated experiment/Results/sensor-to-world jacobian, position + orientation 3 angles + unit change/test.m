min_singular_value_condition_number = [0.1380, 0.1380]
condition_number_condition_number = [0.2022, 0.4894]
configuration = [1, 2]

condition_number_max_singular_value = [2.03e-6 1.09e-6]
condition_number_min_singular_value = [4.10e-7 5.34e-7]
min_singular_value_max_singular_value = [6.19e-6 9.80e-6]
min_singular_value_min_singular_value = [8.56e-7 1.35e-6]

plot(configuration, condition_number_condition_number,'-x', 'LineWidth', 2, 'Color','b')
hold on
plot(configuration, min_singular_value_condition_number,'-x', 'LineWidth', 2, 'Color','r')
title('condition number')
xlim([0 3])
ylim([0 1])
xticks(0:3)
xticklabels({'', '4 3-axis', '12 1-axis', ''})
legend('Maximize rcond', 'Maximize min singular value')
grid on

figure
plot(configuration, condition_number_max_singular_value,'-x', 'LineWidth', 2, 'Color','b')
hold on
plot(configuration, condition_number_min_singular_value,'-x', 'LineWidth', 2, 'Color','c')
hold on 
plot(configuration, min_singular_value_max_singular_value,'-x', 'LineWidth', 2, 'Color','r')
hold on
plot(configuration, min_singular_value_min_singular_value,'-x', 'LineWidth', 2, 'Color','m')

title('max & min singular values')
xlim([0 3])
xticks(0:3)
xticklabels({'', '4 3-axis', '12 1-axis', ''})
legend('Maximize rcond: max singular value', 'Maximize rcond: min singular value', 'Maximize min singular value: max singular value', 'Maximize min singular value: min singular value')

grid on
