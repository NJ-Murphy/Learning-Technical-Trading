function hess = analytical_hess(muest,lamest,sigmaest,thetaest,i,deltanu)
  % covariances
  dmudlam = (-2/sigmaest)*sum(log(i).*(i.^(thetaest-2*lamest)).*(deltanu'-muest*i.^thetaest));
  dsigdlam = -(1/(sigmaest^2))*sum((log(i).*(i.^(-2*lamest))).*(deltanu'-muest*i.^thetaest).^2);
  dsigdmu = -(1/sigmaest^2)*sum((i.^(thetaest-2*lamest)).*(deltanu'-muest*i.^thetaest));
  dthetdlam = ((-2*muest)/sigmaest)*sum((log(i).^2).*(i.^(thetaest-2*lamest)).*(deltanu'-muest*i.^thetaest));
  dmudthet = sum((1/sigmaest)*log(i).*(i.^(thetaest-2*lamest)).*(deltanu'-2*muest*i.^thetaest));
  dthetdsig = (-2*muest/sigmaest^2)*sum(log(i).*(i.^thetaest).*(i.^(-2*lamest)).*(deltanu'-2*muest*i.^thetaest));
  
  % variances 
  d2mu = -(1/sigmaest)*sum(i.^(2*thetaest - 2*lamest));
  d2lam = (-2/sigmaest)*sum((log(i).^2).*(i.^(-2*lamest)).*(deltanu'-muest*i.^thetaest).^2);
  d2thet = (muest/sigmaest)*sum((i.^(-2*lamest)).*(log(i).^2).*(i.^thetaest).*(deltanu'-2*muest*i.^thetaest));
  d2sig = (length(i)/(2*sigmaest^2)) - (1/sigmaest^3)*sum((i.^(-2*lamest)).*(deltanu'-muest*i.^thetaest).^2);
  
  % hessian
  % with theta
  %return(matrix(c(varmu,dmudlam,dmudthet,dsigdmu,dmudlam,varlam,dthetdlam,dsigdlam,dmudthet,dthetdlam,varthet,dthetdsig,dsigdmu,dsigdlam,dthetdsig,varsig),4,4,byrow = T))
  
  %without theta
  hess = vec2mat([d2mu,dmudlam,dsigdmu,dmudlam,d2lam,dsigdlam,dsigdmu,dsigdlam,d2sig],3);
end