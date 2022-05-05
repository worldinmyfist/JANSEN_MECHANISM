function [tpositions, trace] = plot_one_leg(ti, t_estimates, input_angle_step_size, cycles)

close all

t1_0 = t_estimates(1);
t2_0 = t_estimates(2);
t3_0 = t_estimates(3);
t4_0 = t_estimates(4);
t5_0 = t_estimates(5);
t6_0 = t_estimates(6);
t7_0 = t_estimates(7);
t8_0 = t_estimates(8);

trace = zeros(2, floor(2*cycles*pi/input_angle_step_size));

count = 1;
 
for tinput = ti:input_angle_step_size:ti + 2*cycles*pi + input_angle_step_size
    
    [t1, t2, t3, t4, t5, t6, t7, t8] = Jansen_Sequential_Newton_Raphson(tinput, t1_0, t2_0, t3_0, t4_0, t5_0, t6_0, t7_0, t8_0);
    tangles = [tinput, t1, t2, t3, t4, t5, t6, t7, t8];    
    
    tpositions = find_joint_positions(tangles);

    tjointA = [tpositions(1,1) , tpositions(1,2)];
    tjointB = [tpositions(2,1) , tpositions(2,2)];
    tjointC = [tpositions(3,1) , tpositions(3,2)];
    tjointD = [tpositions(4,1) , tpositions(4,2)];
    tjointE = [tpositions(5,1) , tpositions(5,2)];
    tjointF = [tpositions(6,1) , tpositions(6,2)];
    tjoint0 = [tpositions(7,1) , tpositions(7,2)];
    tjoint1 = [tpositions(8,1) , tpositions(8,2)];
    
    drawnow
    clf        
    axis([-100 100 -100 100])
    
    trace(1,count) = tjointF(1);
    trace(2,count) = tjointF(2);
    line(trace(1,1:count), trace(2,1:count));     
    
    patch('Faces', [1 2 3], 'Vertices', [tjoint1; tjointB; tjointC])
    patch('Faces', [1 2 3], 'Vertices', [tjointD; tjointE; tjointF])

    line([tjoint0(1) tjointA(1)], [tjoint0(2) tjointA(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0])
    line([tjointA(1) tjointB(1)], [tjointA(2) tjointB(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjointB(1) tjointC(1)], [tjointB(2) tjointC(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjointE(1) tjointD(1)], [tjointE(2) tjointD(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjointA(1) tjointD(1)], [tjointA(2) tjointD(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjointB(1) tjoint1(1)], [tjointB(2) tjoint1(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjoint1(1) tjointD(1)], [tjoint1(2) tjointD(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjoint1(1) tjointC(1)], [tjoint1(2) tjointC(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30, 'MarkerEdgeColor', [0 0 0])
    line([tjointE(1) tjointF(1)], [tjointE(2) tjointF(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjointD(1) tjointF(1)], [tjointD(2) tjointF(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
    line([tjointC(1) tjointE(1)], [tjointC(2) tjointE(2)], 'Color',[1 0.5 0.5], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)

    t1_0 = t1;
    t2_0 = t2;
    t3_0 = t3;
    t4_0 = t4;
    t5_0 = t5;
    t6_0 = t6;
    t7_0 = t7;
    t8_0 = t8;

    count = count + 1;
    
end
end
