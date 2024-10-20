Try to fix the problem where pred_NL_buckle eigenvalues
don't match the linear eigenvalues.

To debug this - I plan to run nonlinear buckling on the perfect cylinder (no imperfections).
And I should expect the NL eigenvalue predictions to match the linear. Otherwise - I may be
able to figure out the error and correct these predictions.

Check solve vs. solve_local commands to see if the linear eigval is the same (note)
linear eigval w/ solve = 2.502536e+02
linear eigval w/ solve local = 2.9e+07 ? wrong answer..

I fixed that as the linear solve done for the load path was wrong (was just float point nonsense) 
as I didn't solve the stiffness matrix.

Now linear eigval w/ solve local = 500 which is double the value we get with solve.
Try removing the applyBCs on the vars => I suspect this may somewhat mess up the eigvals because
it sets the BCs for lambda = 1 in which may not match out current load point. Could add another optional argument
in with the load factor and to applyBCs on that..

Also removed setBCs in the delta_vars computation (which is probably wrong..) ?
The setBCs calls in to assembler are definitely wrong => but don't seem to significantly affect the eigenvalue
if we already solved the static problem (NL or linear up to that point). Still should fix it just in case..
Try changing the assembleMatType from nonlinearGmat routine.

Another thing for accurate eigvals in nonlinear regime => need to do complex-step on Gmatrix otherwise
the stepsize is 1e-4 and that is way too high for 1e-5 disp control. Need complex-step.

Ok just recovered the correct eigval if I switch back to assembleMat routine (not the assembleNonlinearGmat routine).
Why though? Try printing out delta_vars in the shell element
to see if the Gmatrix computation inputs are different.
Why are all the vars of the load path zero right now? Very confused. Ah => so it was printing out the 0 disp elements first I guess.

The path looks the same from printouts => so is it the initial variables that are different? Why is the eigval 2x for solve_local
with what should be the same inputs? Check the output matrices vs. each other. Oh I just found the mistake, the backwards pert
was using vars[i] not delta_vars[i] (so we had a 2x factor error).

Now I get the same answer lambda = 250, with solve_local as solve for u0 = 0 and du = linear static solve with lambda = 1! And this is using assembleNonlinearGmat too! I might be almost ready to run the nonlinear solve_local in nonlinear buckling again after 
fixing this bug! This error could definitely cause significant errors when vars[i] is very large in the 
nonlinear regime (leading eigvals to be very low).

First, let's do a few sanity checks where if I compute lambda eigval about u0 = 10 * u_LIN and du = u_LIN then I should get eigval -= 10.
Ayy! The pred_NL_buckle works now! Not complete nonsense because of the bug in computing Gmatrix that was using delta_vars and vars pert (really bad bug). Take off the error checking on the eigval probably and stop once the eigenvalue is small enough..
iter, lambda,         |u|/lambda,     dlambda_ds,     loc_eigval,     error,          LIN_buckle,     pred_NL_buckle
11    2.500000e+01    1.496634e-03    1.000000e+00    2.251967e+02    5.950056e-05    2.502536e+02    2.501967e+02 
