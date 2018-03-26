this algorithm works by using fmincon on each element solving for 5 variables (4 lamination parameters and the density) 

in each lambda the lambda is changing so that it cnverges. 

the results showed no convergence and the lambda kept changing without being decreasing 

Dr.Mostafa said this is expected for fmincon since it works by approximating the constraint function as a line this causes the solution to diverge and get out of the bound. 


This is still skeptical since volume fraction and lambda should be decreasing function (this is not the case) 

