function para_final =  gcamp6_model_4para_new(spk, dff, t_frame, para_start)
    
    x.spk           = spk;
    x.dff           = dff;
    x.t_frame       = t_frame;
    mse0            = func_error(x, para_start);
    options         = optimset('Display','off','TolX',1e-5,'TolFun',1e-2);
    para_final      = fminsearch(@(p) func_error(x, p), para_start, options);
end

function mse   = func_error(x, p)
    CaTraces   = func_getCaTraces_general_new(x.spk, x.t_frame, p);
    mse        = sum((CaTraces - x.dff).^2);
end