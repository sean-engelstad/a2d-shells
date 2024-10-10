From TACSShellElement.h:
    from addResidual() method:
    
    // can't do strain energy with A2D unless you put A2D inside ShellElementModel also
    // // strain energy stack (backprop to disp grads and tying strain later)
    // auto strain_energy_stack = A2D::MakeStack(
    //   // intermediate steps for u0x, u1x, etc. on TODO
    //   // last step
    //   A2D::SymMatMultTrace(E, S, ES_prod),
    //   A2D::Eval(0.5 * detXd * ES_prod, Uelem)
    // );