//[[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>
using namespace arma;

//[[Rcpp::export()]]
arma::vec S1_wls(arma::vec summary, double n){
  double w1=4/(4+pow(n, 0.75));
  double w3=1;
  arma::vec beta_hat(2, arma::fill::zeros);
  
  double mu_hat=w1*(summary(0)+summary(2))/2+(1-w1)*summary(1);
  double sigma_hat=w3*(summary(2)-summary(0))/(2*R::qnorm((n-0.375)/(n+0.25), 0, 1, 1, 0));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S2_wls(arma::vec summary, double n){
  double w2=0.7+0.39/n;
  double w3=0;
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w2*(summary(0)+summary(2))/2+(1-w2)*summary(1);
  double sigma_hat=(1-w3)*(summary(2)-summary(0))/(2*R::qnorm((0.75*n-0.125)/(n+0.25), 0, 1, 1, 0));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S3_wls(arma::vec summary, double n){
  double w1=2.2/(2.2+pow(n, 0.75));
  double w2=0.7-0.72/pow(n, 0.55);
  double w3=1/(1+0.07*pow(n, 0.6));
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w1*(summary(0)+summary(4))/2+w2*(summary(1)+summary(3))/2+(1-w1-w2)*summary(2);
  double sigma_hat=w3*(summary(4)-summary(0))/(2*R::qnorm((n-0.375)/(n+0.25), 0, 1, 1, 0))+(1-w3)*(summary(3)-summary(1))/(2*R::qnorm((0.75*n-0.125)/(n+0.25), 0, 1, 1, 0));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

double Q1(double r, int n){
  double gamma=-R::qnorm(r/(n+1), 0, 1, 1, 0)/sqrt(2);
  return sqrt(2*datum::pi)*exp(gamma*gamma);
}

double Q2(double r, int n){
  double gamma=-R::qnorm(r/(n+1), 0, 1, 1, 0)/sqrt(2);
  return -2*sqrt(2)*datum::pi*gamma*exp(2*gamma*gamma);
}

double Q3(double r, int n){
  double gamma=-R::qnorm(r/(n+1), 0, 1, 1, 0)/sqrt(2);
  return 2*sqrt(2)*pow(datum::pi, 1.5)*exp(3*gamma*gamma)*(1+4*gamma*gamma);
}

double Q4(double r, int n){
  double gamma=-R::qnorm(r/(n+1), 0, 1, 1, 0)/sqrt(2);
  return -4*sqrt(2)*pow(datum::pi, 2)*gamma*exp(4*gamma*gamma)*(7+12*gamma*gamma);
}

double Expectation_yang(double r, int n){
  double pr=r/(n+1);
  double qr=1-pr;
  double value=(qr-pr)*Q3(r, n)/3+pr*qr*Q4(r, n)/8;
  
  return R::qnorm(pr, 0, 1, 1, 0)+pr*qr*Q2(r, n)/(2*n+4)+pr*qr*value/((n+2)*(n+2));
}

double Covariance_yang(double r, double s, int n){
  double pr=r/(n+1);
  double qr=1-pr;
  double ps=s/(n+1);
  double qs=1-ps;
  double value1=(qr-pr)*Q2(r, n)*Q1(s, n)+(qs-ps)*Q1(r, n)*Q2(s, n);
  double value2=pr*qr*Q3(r, n)*Q1(s, n)/2+ps*qs*Q3(s, n)*Q1(r, n)/2+pr*qs*Q2(r, n)*Q2(s, n)/2;
  
  return pr*qs*Q1(r, n)*Q1(s, n)/(n+2)+pr*qs*(value1+value2)/((n+2)*(n+2));
}

arma::mat alpha_yang(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat alpha=ones(2, len);
  for(int i=0; i<len; i++){
    alpha(1, i)=Expectation_yang(index(i), n);
  }
  return alpha.t();
}

arma::mat cov_yang(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat cov=zeros(len, len);
  for(int i=0; i<len; i++){
    for(int j=0; j<len; j++)
      if(i>j){cov(i, j)=cov(j, i);}
      else{cov(i, j)=Covariance_yang(index(i), index(j), n);}
  }
  return cov;
}

//[[Rcpp::export()]]
arma::vec S1_yang(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/2)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang(index, n);
  cov=cov_yang(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S2_yang(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang(index, n);
  cov=cov_yang(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S3_yang(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang(index, n);
  cov=cov_yang(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec tertiles_yang(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/3)+1, floor(2*n/3)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang(index, n);
  cov=cov_yang(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec quintiles_yang(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/5)+1, floor(2*n/5)+1, floor(3*n/5)+1, floor(4*n/5)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang(index, n);
  cov=cov_yang(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec deciles_yang(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/10)+1, floor(2*n/10)+1, floor(3*n/10)+1, floor(4*n/10)+1, floor(5*n/10)+1, floor(6*n/10)+1, floor(7*n/10)+1, floor(8*n/10)+1, floor(9*n/10)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang(index, n);
  cov=cov_yang(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

double Q1_logit(double r, int n){
  double pr=r/(n+1);
  double s=sqrt(3)/datum::pi;
  return s/(pr*(1-pr));
}

double Q2_logit(double r, int n){
  double pr=r/(n+1);
  double s=sqrt(3)/datum::pi;
  return s*(2*pr-1)/(pr*pr*(1-pr)*(1-pr));
}

double Q3_logit(double r, int n){
  double pr=r/(n+1);
  double s=sqrt(3)/datum::pi;
  return 2*s*(3*pr*pr-3*pr+1)/(pr*pr*pr*(1-pr)*(1-pr)*(1-pr));
}

double Q4_logit(double r, int n){
  double pr=r/(n+1);
  double s=sqrt(3)/datum::pi;
  return 6*s*(4*pr*pr*pr-6*pr*pr+4*pr-1)/(pr*pr*pr*pr*(1-pr)*(1-pr)*(1-pr)*(1-pr));
}

double Expectation_yang_logit(double r, int n){
  double pr=r/(n+1);
  double qr=1-pr;
  double s=sqrt(3)/datum::pi;
  double value=(qr-pr)*Q3_logit(r, n)/3+pr*qr*Q4_logit(r, n)/8;
  
  return 0+s*log(pr/qr)+pr*qr*Q2_logit(r, n)/(2*n+4)+pr*qr*value/((n+2)*(n+2));
}

double Covariance_yang_logit(double r, double s, int n){
  double pr=r/(n+1);
  double qr=1-pr;
  double ps=s/(n+1);
  double qs=1-ps;
  double value1=(qr-pr)*Q2_logit(r, n)*Q1_logit(s, n)+(qs-ps)*Q1_logit(r, n)*Q2_logit(s, n);
  double value2=pr*qr*Q3_logit(r, n)*Q1_logit(s, n)/2+ps*qs*Q3_logit(s, n)*Q1_logit(r, n)/2+pr*qs*Q2_logit(r, n)*Q2_logit(s, n)/2;
  
  return pr*qs*Q1_logit(r, n)*Q1_logit(s, n)/(n+2)+pr*qs*(value1+value2)/((n+2)*(n+2));
}

arma::mat alpha_yang_logit(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat alpha=ones(2, len);
  for(int i=0; i<len; i++){
    alpha(1, i)=Expectation_yang_logit(index(i), n);
  }
  return alpha.t();
}

arma::mat cov_yang_logit(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat cov=zeros(len, len);
  for(int i=0; i<len; i++){
    for(int j=0; j<len; j++)
      if(i>j){cov(i, j)=cov(j, i);}
      else{cov(i, j)=Covariance_yang_logit(index(i), index(j), n);}
  }
  return cov;
}

//[[Rcpp::export()]]
arma::vec S1_yang_logit(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/2)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_logit(index, n);
  cov=cov_yang_logit(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S2_yang_logit(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_logit(index, n);
  cov=cov_yang_logit(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S3_yang_logit(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_logit(index, n);
  cov=cov_yang_logit(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec tertiles_yang_logit(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/3)+1, floor(2*n/3)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_logit(index, n);
  cov=cov_yang_logit(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec quintiles_yang_logit(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/5)+1, floor(2*n/5)+1, floor(3*n/5)+1, floor(4*n/5)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_logit(index, n);
  cov=cov_yang_logit(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec deciles_yang_logit(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/10)+1, floor(2*n/10)+1, floor(3*n/10)+1, floor(4*n/10)+1, floor(5*n/10)+1, floor(6*n/10)+1, floor(7*n/10)+1, floor(8*n/10)+1, floor(9*n/10)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_logit(index, n);
  cov=cov_yang_logit(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

double Q1_lap(double r, int n){
  double pr=r/(n+1);
  if(pr>=0.5){return 1/(sqrt(2)*(1-pr));}
  else{return 1/(sqrt(2)*pr);}
}

double Q2_lap(double r, int n){
  double pr=r/(n+1);
  if(pr>=0.5){return 1/(sqrt(2)*(1-pr)*(1-pr));}
  else{return -1/(sqrt(2)*pr*pr);}
}

double Q3_lap(double r, int n){
  double pr=r/(n+1);
  if(pr>=0.5){return sqrt(2)/((1-pr)*(1-pr)*(1-pr));}
  else{return sqrt(2)/(pr*pr*pr);}
}

double Q4_lap(double r, int n){
  double pr=r/(n+1);
  if(pr>=0.5){return 3*sqrt(2)/((1-pr)*(1-pr)*(1-pr)*(1-pr));}
  else{return -3*sqrt(2)/(pr*pr*pr*pr);}
}

double Expectation_yang_lap(double r, int n){
  double pr=r/(n+1);
  double qr=1-pr;
  double value=(qr-pr)*Q3_lap(r, n)/3+pr*qr*Q4_lap(r, n)/8;
  
  return 0-sign(pr-0.5)*log(1-2*abs(pr-0.5))/sqrt(2)+pr*qr*Q2_lap(r, n)/(2*n+4)+pr*qr*value/((n+2)*(n+2));
}

double Covariance_yang_lap(double r, double s, int n){
  double pr=r/(n+1);
  double qr=1-pr;
  double ps=s/(n+1);
  double qs=1-ps;
  double value1=(qr-pr)*Q2_lap(r, n)*Q1_lap(s, n)+(qs-ps)*Q1_lap(r, n)*Q2_lap(s, n);
  double value2=pr*qr*Q3_lap(r, n)*Q1_lap(s, n)/2+ps*qs*Q3_lap(s, n)*Q1_lap(r, n)/2+pr*qs*Q2_lap(r, n)*Q2_lap(s, n)/2;
  
  return pr*qs*Q1_lap(r, n)*Q1_lap(s, n)/(n+2)+pr*qs*(value1+value2)/((n+2)*(n+2));
}

arma::mat alpha_yang_lap(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat alpha=ones(2, len);
  for(int i=0; i<len; i++){
    alpha(1, i)=Expectation_yang_lap(index(i), n);
  }
  return alpha.t();
}

arma::mat cov_yang_lap(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat cov=zeros(len, len);
  for(int i=0; i<len; i++){
    for(int j=0; j<len; j++)
      if(i>j){cov(i, j)=cov(j, i);}
      else{cov(i, j)=Covariance_yang_lap(index(i), index(j), n);}
  }
  return cov;
}

//[[Rcpp::export()]]
arma::vec S1_yang_lap(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/2)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_lap(index, n);
  cov=cov_yang_lap(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S2_yang_lap(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_lap(index, n);
  cov=cov_yang_lap(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S3_yang_lap(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_lap(index, n);
  cov=cov_yang_lap(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec tertiles_yang_lap(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/3)+1, floor(2*n/3)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_lap(index, n);
  cov=cov_yang_lap(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec quintiles_yang_lap(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/5)+1, floor(2*n/5)+1, floor(3*n/5)+1, floor(4*n/5)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_lap(index, n);
  cov=cov_yang_lap(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec deciles_yang_lap(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/10)+1, floor(2*n/10)+1, floor(3*n/10)+1, floor(4*n/10)+1, floor(5*n/10)+1, floor(6*n/10)+1, floor(7*n/10)+1, floor(8*n/10)+1, floor(9*n/10)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_yang_lap(index, n);
  cov=cov_yang_lap(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

double Expectation_bala(double r, int n){
  double pr=r/(n+1);
  
  return R::qnorm(pr, 0, 1, 1, 0)+Q2(r, n)*pr*((r+1)/(n+2)-pr)/2;
}


double Covariance_bala(double r, double s, int n){
  if(((r==1)&&(s==1))||((r==n)&&(s==n))){
    return pow(-0.1565*log(log(n))+0.8949, 4);
  }
  else{
    // the following codes may contain errors in the original paper 
    // by Balakrishnan et al. (2022)
    // so we replaced them by his taylor expansion as a compromise
    // double q1=r/n;
    // double q2=s/n;
    // double a=sqrt(1-4*(q1-0.5)*(q2-0.5)+2*pow(abs(q1-q2), 0.5))/8-0.51;
    // double b=sqrt(1-4*(q1-0.5)*(q2-0.5)+2*pow(abs(q1-q2), 0.5))/2.8+1.3;
    // 
    // return pow(a*log(log(n))+b, 4);
    
    double pr=r/(n+1);
    double qs=s/(n+1);
    double term1=R::qnorm(pr, 0, 1, 1, 0)*R::qnorm(qs, 0, 1, 1, 0);
    double term2=Q2(r, n)*R::qnorm(qs, 0, 1, 1, 0)*pr*((r+1)/(n+2)-pr)/2;
    double term3=Q1(r, n)*Q1(s, n)*pr*((s+1)/(n+2)-qs);
    double term4=R::qnorm(pr, 0, 1, 1, 0)*Q2(s, n)*qs*((s+1)/(n+2)-qs)/2;
    
    return term1+term2+term3+term4-Expectation_bala(r, n)*Expectation_bala(s, n);
    
  }
}

arma::mat alpha_bala(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat alpha=ones(2, len);
  for(int i=0; i<len; i++){
    alpha(1, i)=Expectation_bala(index(i), n);
  }
  return alpha.t();
}

arma::mat cov_bala(arma::vec index, int n){
  int len=index.n_elem;
  arma::mat cov=zeros(len, len);
  for(int i=0; i<len; i++){
    for(int j=0; j<len; j++)
      if(i>j){cov(i, j)=cov(j, i);}
      else{cov(i, j)=Covariance_bala(index(i), index(j), n);}
  }
  return cov;
}

//[[Rcpp::export()]]
arma::vec S1_bala(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/2)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_bala(index, n);
  cov=cov_bala(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S2_bala(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_bala(index, n);
  cov=cov_bala(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}


//[[Rcpp::export()]]
arma::vec S3_bala(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={1, floor(n/4)+1, floor(n/2)+1, floor(3*n/4)+1, n};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_bala(index, n);
  cov=cov_bala(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec tertiles_bala(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/3)+1, floor(2*n/3)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_bala(index, n);
  cov=cov_bala(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec quintiles_bala(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/5)+1, floor(2*n/5)+1, floor(3*n/5)+1, floor(4*n/5)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_bala(index, n);
  cov=cov_bala(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec deciles_bala(arma::vec summary, double n){
  int len=summary.n_elem;
  arma::vec index={floor(n/10)+1, floor(2*n/10)+1, floor(3*n/10)+1, floor(4*n/10)+1, floor(5*n/10)+1, floor(6*n/10)+1, floor(7*n/10)+1, floor(8*n/10)+1, floor(9*n/10)+1};
  arma::mat cov=zeros(len, len);
  arma::mat alpha=zeros(2, len);
  
  alpha=alpha_bala(index, n);
  cov=cov_bala(index, n);
  
  arma::mat omega=inv(cov);
  
  return inv(alpha.t()*omega*alpha)*alpha.t()*omega*summary;
}

//[[Rcpp::export()]]
arma::vec S1_lap(arma::vec summary, double n){
  double w1=1/(-0.468533+0.598561*pow(n, 1.60567));
  double w3=1;
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w1*(summary(0)+summary(2))/2+(1-w1)*summary(1);
  double sigma_hat=w3*(summary(2)-summary(0))/(2/(-0.0084688+1.46221*pow(log(n), -0.988499)));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S2_lap(arma::vec summary, double n){
  double w2=1/(-0.0449944+0.870993*pow(n, 0.555273));
  double w3=0;
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w2*(summary(0)+summary(2))/2+(1-w2)*summary(1);
  double sigma_hat=(1-w3)*(summary(2)-summary(0))/(2/(2.0424+2.91911*pow(n, -1.20853)));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S3_lap(arma::vec summary, double n){
  double w1=1.12159*exp(-1.13711-0.476207*n);
  double w2=1/(0.297173+0.759076*pow(n, 0.576925));
  double w3=9.30594/(10.8703+0.064734*n+0.581161*pow(log(n), 1.95934));
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w1*(summary(0)+summary(4))/2+w2*(summary(1)+summary(3))/2+(1-w1-w2)*summary(2);
  double sigma_hat=w3*(summary(4)-summary(0))/(2/(-0.0084688+1.46221*pow(log(n), -0.988499)))+(1-w3)*(summary(3)-summary(1))/(2/(2.0424+2.91911*pow(n, -1.20853)));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S1_logit(arma::vec summary, double n){
  double w1=65.4517/(50.3177+27.9285*pow(n, 0.992862));
  double w3=1;
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w1*(summary(0)+summary(2))/2+(1-w1)*summary(1);
  double sigma_hat=w3*(summary(2)-summary(0))/(2/(0.182894+1.19771*pow(n, -0.426478)));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S2_logit(arma::vec summary, double n){
  double w2=1/(1.66681-1.05697*pow(n, -0.996803));
  double w3=0;
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w2*(summary(0)+summary(2))/2+(1-w2)*summary(1);
  double sigma_hat=(1-w3)*(summary(2)-summary(0))/(2/(1.65277+3.1854*pow(n, -1.13223)));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}

//[[Rcpp::export()]]
arma::vec S3_logit(arma::vec summary, double n){
  double w1=109.931/(-17.958+117.682*pow(n, 0.992282));
  double w2=(0.604702*n)/(0.991393+pow(n, 1.00126));
  double w3=448.738/((404.359+16.8838*pow(n, 0.660158))*pow(n, 0.143355));
  arma::vec beta_hat(2, fill::zeros);
  
  double mu_hat=w1*(summary(0)+summary(4))/2+w2*(summary(1)+summary(3))/2+(1-w1-w2)*summary(2);
  double sigma_hat=w3*(summary(4)-summary(0))/(2/(0.182894+1.19771*pow(n, -0.426478)))+(1-w3)*(summary(3)-summary(1))/(2/(1.65277+3.1854*pow(n, -1.13223)));
  
  beta_hat(0)=mu_hat;
  beta_hat(1)=sigma_hat;
  
  return beta_hat;
}



