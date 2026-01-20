function Res = Res_vp ( dg_ , F__ , Gve_ , K__ , sin_phi_ , sin_psi , eta_vp , dt) 
Res = 0 - F__ + dg_*Gve_ - dg_*K__.*sin_phi_.*sin_psi - dg_*eta_vp/dt;
end