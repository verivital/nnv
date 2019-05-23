% execute all zonotope tests
t_pause = 0.1;


close all ; clearAllExceptVars(t_pause);
test_zono_convexHull
pause(t_pause);

close all ; clearAllExceptVars(t_pause);
test_zono_convexHull_with_linearTransform
pause(t_pause);

close all ; clearAllExceptVars(t_pause);
test_zono_getOrientedBox
pause(t_pause);

close all ; clearAllExceptVars(t_pause);
test_zono_getSupInfinityNorm
pause(t_pause);

close all ; clearAllExceptVars(t_pause);
test_zono_getVertices
pause(t_pause);

close all ; clearAllExceptVars(t_pause);
test_zono_orderReduction
pause(t_pause);

close all ; clearAllExceptVars(t_pause);
zono_test
pause(t_pause);

close all ; clear all;