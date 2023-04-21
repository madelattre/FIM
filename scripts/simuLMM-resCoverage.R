

## expected coverage rate
rate <- 0.95 


# ## Print the coverage rates 
# ## ## Coverage rates of confidence ellipsoids ...
# theta.est <- matrix(c(est.mle@beta,variances[1,'vcov'],variances[2,'vcov']),ncol=1)
# unit.vect <- matrix(rep(1,length(theta.est)),ncol=1)
# 
# ## ... based on the true Fisher information matrix
# 
# #conf.inf.fisher.true <- theta.est -
# #  1.96/sqrt(n)*sqrt(diag(solve(fisher)))
# #conf.sup.fisher.true <- theta.est +
# #  1.96/sqrt(n)*sqrt(diag(solve(fisher)))
# 
# #coverage.true.fisher <- coverage.true.fisher +
# #  (theta.true<=conf.sup.fisher.true)*(theta.true>=conf.inf.fisher.true)
# coverage.true.fisher <- coverage.true.fisher + 
#   (n*t(as.matrix(theta.true-theta.est))%*%fisher%*%as.matrix(theta.true-theta.est)
#    <= qchisq(rate,3))
# 
# ## ... based on Isco computed in the MLE value of the parameters
# # conf.inf.isco.theta.est <- theta.est -
# #   1.96/sqrt(n)*sqrt(diag(solve(resIsco.theta.est[,,k])))
# # conf.sup.isco.theta.est <- theta.est +
# #   1.96/sqrt(n)*sqrt(diag(solve(resIsco.theta.est[,,k])))
# # 
# # coverage.isco.theta.est <- coverage.isco.theta.est +
# #   (theta.true<=conf.sup.isco.theta.est)*(theta.true>=conf.inf.isco.theta.est)
# coverage.isco.theta.est <- coverage.isco.theta.est +
#   (n*t(as.matrix(theta.true-theta.est))%*%resIsco.theta.est[,,k]%*%as.matrix(theta.true-theta.est)
#    <= qchisq(rate,3))
# 
# ## ... based on Isco computed in the true parameter values
# # conf.inf.isco.theta.true <- theta.est -
# #   1.96/sqrt(n)*sqrt(diag(solve(resIsco.theta.true[,,k])))
# # conf.sup.isco.theta.true <- theta.est +
# #   1.96/sqrt(n)*sqrt(diag(solve(resIsco.theta.true[,,k])))
# # 
# # coverage.isco.theta.true <- coverage.isco.theta.true +
# #   (theta.true<=conf.sup.isco.theta.true)*(theta.true>=conf.inf.isco.theta.true)
# coverage.isco.theta.true <- coverage.isco.theta.true +
#   (n*t(as.matrix(theta.true-theta.est))%*%resIsco.theta.true[,,k]%*%as.matrix(theta.true-theta.est)
#    <= qchisq(rate,3))
# 
# ## ... based on Iobs computed in the MLE value of the parameters
# # conf.inf.iobs.theta.est <- theta.est -
# #   1.96/sqrt(n)*sqrt(diag(solve(resIobs.theta.est[,,k])))
# # conf.sup.iobs.theta.est <- theta.est +
# #   1.96/sqrt(n)*sqrt(diag(solve(resIobs.theta.est[,,k])))
# # 
# # coverage.iobs.theta.est <- coverage.iobs.theta.est +
# #   (theta.true<=conf.sup.iobs.theta.est)*(theta.true>=conf.inf.iobs.theta.est)
# coverage.iobs.theta.est <- coverage.iobs.theta.est +
#   (n*t(as.matrix(theta.true-theta.est))%*%resIobs.theta.est[,,k]%*%as.matrix(theta.true-theta.est)
#    <= qchisq(rate,3))
# 
# ## ... based on Iobs computed in the true parameter values
# # conf.inf.iobs.theta.true <- theta.est -
# #   1.96/sqrt(n)*sqrt(diag(solve(resIobs.theta.true[,,k])))
# # conf.sup.iobs.theta.true <- theta.est +
# #   1.96/sqrt(n)*sqrt(diag(solve(resIobs.theta.true[,,k])))
# # 
# # coverage.iobs.theta.true <- coverage.iobs.theta.true +
# #   (theta.true<=conf.sup.iobs.theta.true)*(theta.true>=conf.inf.iobs.theta.true)
# coverage.iobs.theta.true <- coverage.iobs.theta.true +
#   (n*t(as.matrix(theta.true-theta.est))%*%resIobs.theta.true[,,k]%*%as.matrix(theta.true-theta.est)
#    <= qchisq(rate,3))  
# ## MD, cat ou print de ce qu'on conserve
# print(coverage.isco.theta.true/nsim*100)
# # print(coverage.isco.theta.est/nsim*100)
# # print(coverage.iobs.theta.true/nsim*100)
# # print(coverage.iobs.theta.est/nsim*100)
# # print(coverage.true.fisher/nsim*100)
# # 