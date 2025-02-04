clc
clear all
close all

mA=1;mB=1;mC=1;mD=1;l1=1;l2=1;l3=1;l4=1;l0=1;

syms d2q1 d2q2 d2q3 d2q4   
syms dq1 dq2 dq3 dq4
syms q1 q2 q3 q4
syms la1 la2
syms T1 T2


Term1=(mA/3+mB)*l1^2*d2q1 + mB*l1*l2/4 * d2q2 *cos(q1-q2) - mB*l1^l2/4*(dq1-dq2)*dq2*sin(q1-q2) + ...
    +mB*l1*l2/4*dq1*dq2*sin(q1-q2)-la1*l1*sin(q1)+la2*l1*cos(q1)-T1

Term2=mB*l2^2/3*d2q2 + mB*l1*l2/4*d2q1*cos(q1-q2)-mB*l1*l2/4*(dq1-dq2)*dq1*sin(q1-q2)-mB*l1*l2/4*dq1*dq2*sin(q1-q2) ...
    -la1*l2*sin(q2)+l2*la2*cos(q2)

Term3=mC*l3^2/3*d2q3+mC*l3*l4/4*d2q4*cos(q3-q4)-mC*l3*l4/4*(dq3-dq4)*dq4*sin(q3-q4)-mC*l3*l4/4*dq3*dq4*sin(q3-q4)+...
    la1*l3*sin(q3)-la2*l3*cos(q3)

Term4=(mD/3+mC)*l4^2*d2q4+mC*l3*l4/4*d2q3*cos(q3-q4)-mC*l3*l4/4*(dq3-dq4)*dq3*sin(q3-q4)+mC*l3*l4/4*dq3*dq4*sin(q3-q4) ...
    +la1*l4*sin(q4)-l4*la2*cos(q4)-T2                                                         %%Equations of motion derived from Lagrangian equation plus lagrange multipliers

Term5=l1*cos(q1)+l2*cos(q2)-l3*cos(q3)-l4*cos(q4)-l0                                          %%Constraint 1

Term6=l1*sin(q1)+l2*sin(q2)-l3*sin(q3)-l4*sin(q4)                                             %%Constraint 2

dTerm5=-l1*dq1*sin(q1)-l2*dq2*sin(q2)+l3*dq3*sin(q3)+l4*dq4*sin(q4)

dTerm6=l1*dq1*cos(q1)+l2*dq2*cos(q2)-l3*dq3*cos(q3)-l4*dq4*cos(q4)                            %%First derivatives of constraints

d2Term5=-l1*d2q1*sin(q1)-l1*dq1^2*cos(q1)-l2*d2q2*sin(q2)-l2*dq2^2*cos(q2)...
    +l3*d2q3*sin(q3)+l3*dq3^2*cos(q3)+l4*d2q4*sin(q4)+l4*dq4^2*cos(q4)

d2Term6=l1*d2q1*cos(q1)-l1*dq1^2*sin(q1)+l2*d2q2*cos(q2)-l2*dq2^2*sin(q2)-l3*d2q3*cos(q3)...  %%Second
    +l3*dq3^2*sin(q3)-l4*d2q4*cos(q4)+l4*dq4^2*sin(q4)

solve(Term1,Term2,Term3,Term4,d2Term5,d2Term6,d2q1,d2q2,d2q3,d2q4,la1,la2)                    %%With four lagrangian equation and second derivatives of constraints, we can solve for d2q1-d2q4

