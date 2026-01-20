M11_ = 1 + dg_.*(D11.*d2Qdsxxdsxx_+D12.*d2Qdsyydsxx_+D13.*d2Qdsxydsxx_+D14.*d2Qdszzdsxx_);
M12_ = 0 + dg_.*(D11.*d2Qdsxxdsyy_+D12.*d2Qdsyydsyy_+D13.*d2Qdsxydsyy_+D14.*d2Qdszzdsyy_);
M13_ = 0 + dg_.*(D11.*d2Qdsxxdsxy_+D12.*d2Qdsyydsxy_+D13.*d2Qdsxydsxy_+D14.*d2Qdszzdsxy_);
M14_ = 0 + dg_.*(D11.*d2Qdsxxdszz_+D12.*d2Qdsyydszz_+D13.*d2Qdsxydszz_+D14.*d2Qdszzdszz_);

M21_ = 0 + dg_.*(D21.*d2Qdsxxdsxx_+D22.*d2Qdsyydsxx_+D23.*d2Qdsxydsxx_+D24.*d2Qdszzdsxx_);
M22_ = 1 + dg_.*(D21.*d2Qdsxxdsyy_+D22.*d2Qdsyydsyy_+D23.*d2Qdsxydsyy_+D24.*d2Qdszzdsyy_);
M23_ = 0 + dg_.*(D21.*d2Qdsxxdsxy_+D22.*d2Qdsyydsxy_+D23.*d2Qdsxydsxy_+D24.*d2Qdszzdsxy_);
M24_ = 0 + dg_.*(D21.*d2Qdsxxdszz_+D22.*d2Qdsyydszz_+D23.*d2Qdsxydszz_+D24.*d2Qdszzdszz_);

M31_ = 0 + dg_.*(D31.*d2Qdsxxdsxx_+D32.*d2Qdsyydsxx_+D33.*d2Qdsxydsxx_+D34.*d2Qdszzdsxx_);
M32_ = 0 + dg_.*(D31.*d2Qdsxxdsyy_+D32.*d2Qdsyydsyy_+D33.*d2Qdsxydsyy_+D34.*d2Qdszzdsyy_);
M33_ = 1 + dg_.*(D31.*d2Qdsxxdsxy_+D32.*d2Qdsyydsxy_+D33.*d2Qdsxydsxy_+D34.*d2Qdszzdsxy_);
M34_ = 0 + dg_.*(D31.*d2Qdsxxdszz_+D32.*d2Qdsyydszz_+D33.*d2Qdsxydszz_+D34.*d2Qdszzdszz_);

M41_ = 0 + dg_.*(D41.*d2Qdsxxdsxx_+D42.*d2Qdsyydsxx_+D43.*d2Qdsxydsxx_+D44.*d2Qdszzdsxx_);
M42_ = 0 + dg_.*(D41.*d2Qdsxxdsyy_+D42.*d2Qdsyydsyy_+D43.*d2Qdsxydsyy_+D44.*d2Qdszzdsyy_);
M43_ = 0 + dg_.*(D41.*d2Qdsxxdsxy_+D42.*d2Qdsyydsxy_+D43.*d2Qdsxydsxy_+D44.*d2Qdszzdsxy_);
M44_ = 1 + dg_.*(D41.*d2Qdsxxdszz_+D42.*d2Qdsyydszz_+D43.*d2Qdsxydszz_+D44.*d2Qdszzdszz_);


[M_inv1_1,M_inv1_2,M_inv1_3,M_inv1_4,M_inv2_1,M_inv2_2,M_inv2_3,M_inv2_4,M_inv3_1,M_inv3_2,M_inv3_3,M_inv3_4,M_inv4_1,M_inv4_2,M_inv4_3,M_inv4_4] = ...
get_inv_4x4 ( M11_ , M21_  , M31_ , M41_ , M12_ , M22_  , M32_  , M42_ , M13_ , M23_  , M33_  , M43_ , M14_ , M24_  , M34_  , M44_ ) ;
M_inv_ = [M_inv1_1 M_inv1_2 M_inv1_3 M_inv1_4 M_inv2_1 M_inv2_2 M_inv2_3 M_inv2_4 M_inv3_1 M_inv3_2 M_inv3_3 M_inv3_4 M_inv4_1 M_inv4_2 M_inv4_3 M_inv4_4];
