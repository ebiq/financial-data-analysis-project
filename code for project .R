data("Nile")
library(dlm)
data_new <- window(Nile, start = 1880, end = 1950)
buildNile <- function(x) {
  dlmModPoly(order = 1, dV = x[1], dW = x[2])
}
mleNile <- dlmMLE(
  y = data_new,
  parm = rep(100, 2),
  build = buildNile,
  lower = rep(1e-8, 2)
)
cat("Estimated σ²_ε:", round(mleNile$par[1], 4), "\n")
cat("Estimated σ²_η:", round(mleNile$par[2], 4), "\n")
cat("Log-likelihood value:", round(-mleNile$value, 4), "\n")
modNile <- dlmModPoly(order = 1, dV = mleNile$par[1],
                      dW = mleNile$par[2], m0 = 0, C0 = 1e8)
filtered <- dlmFilter(data_new, modNile)
filtered$se <- sqrt(unlist(dlmSvd2var(filtered$U.C, filtered$D.C)))







plot(data_new, lty = 2, ylab = "Flow", xlab = "Year", main = "Nile River Flow with 95% Confidence Interval")  
lines(dropFirst(filtered$m), col = "blue", lwd = 2)  
lines(dropFirst(filtered$m - 1.96 * filtered$se), lty = 3, col = "red")  
lines(dropFirst(filtered$m + 1.96 * filtered$se), lty = 3, col = "red")  
points(time(data_new), data_new, pch = 16, col = "black", cex = 0.6)  




f <- filtered$f      
y <- filtered$y       
  
vt <- y - f
plot(vt, type = "l", col = "black", lwd = 1.5,
     main = "Prediction Errors (vt)", xlab = "Year", ylab = "Error",
     ylim = c(-500, 500))  
abline(h = 0, col = "red", lty = 2)








smoothed <- dlmSmooth(filtered)


plot(data_new, type = "o", col = "black", pch = 16, cex = 0.6,
     xlab = "Year", ylab = "Flow", main = "Nile River Flow: Original, Filtered, and Smoothed Data")  

lines(dropFirst(filtered$m), col = "blue", lwd = 2, lty = 1)  
lines(dropFirst(smoothed$s), col = "red", lwd = 2, lty = 2)  






























#glgssm




  
  
  data("Seatbelts")

  selected_data <- window(Seatbelts[, "DriversKilled"], start=c(1974,1), end=c(1984,12))
  
  log_KSI <- log(selected_data)
  ts.plot(log_KSI, col="black", lty=1, xlab="Year", ylab="Log(Deaths)",
          main="Driver killed in UK Road Accidents (1974-1984)")
  points(time(log_KSI), log_KSI, pch=16, cex=0.6) 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  drivers_ts <- ts(Seatbelts[, "drivers"], start = c(1969, 1), frequency = 12)
  data_new <- log(window(drivers_ts, start = c(1974, 1), end = c(1984, 12)))
  
  dlmModel <- dlmModPoly(order = 1) + dlmModSeas(frequency = 12)

  buildFun <- function(x) {
    diag(W(dlmModel))[1] <- exp(x[1])  
    diag(W(dlmModel))[2:12] <- exp(x[2])  
    V(dlmModel) <- exp(x[3])  
    return(dlmModel)
  }
  
 
  fit <- dlmMLE(data_new, parm = rep(0, 3), build = buildFun, method = "L-BFGS-B")
  
 
  dlmFinal <- buildFun(fit$par)
  

  sigma_1 <- drop(V(dlmFinal))  
  sigma_2 <- diag(W(dlmFinal))[1]  
  sigma_3 <- diag(W(dlmFinal))[2]  
  

  variance_results <- data.frame(
    Component = c("Local Level", "Seasonal", "Irregular"),
    Value = c(sigma_2, sigma_3, sigma_1)
  )
  
  print(variance_results)
  
 
  smoothed <- dlmSmooth(data_new, mod = dlmFinal)
  

  trend <- dropFirst(smoothed$s[, 1]) 
  seasonal <- dropFirst(smoothed$s[, 2])  
  irregular <- data_new - trend - seasonal  
  

  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  
  # (i) Original Data + Smoothed Trend
  plot(data_new, type = "l", col = "black", ylab = "Log(Drivers)", main = "(i) Original Data + Smoothed Trend")
  lines(trend, col = "blue", lwd = 2)
  
  # (ii) Seasonal Component
  plot(seasonal, type = "l", col = "black", ylab = "Seasonal", main = "(ii) Seasonal Component")
  abline(h = 0, col = "blue")  
  # (iii) Irregular Component
  plot(irregular, type = "l", col = "black", ylab = "Irregular", main = "(iii) Irregular Component")
  abline(h = 0, col = "blue") 
  
  par(mfrow = c(1, 1))  
  
  
  
  
  
  
  
  
  library(datasets)
  data("Seatbelts") 
  drivers_ts <- ts(Seatbelts[, "drivers"], start = c(1969, 1), frequency = 12)
  

  drivers_window <- window(drivers_ts, start = c(1974, 1), end = c(1984, 12))
  

  log_drivers <- log(drivers_window)
  

  decomp_log_drivers <- decompose(log_drivers, type = "additive")
    

    par(mfrow = c(3, 1))  
    
   
    plot(log_drivers, type = "l", col = "black", ylab = "Log(Drivers)", main = "(i) Original Data +  Trend")
    lines(decomp_log_drivers$trend, col = "blue", lwd = 2)
    

    plot(decomp_log_drivers$seasonal, type = "l", col = "black", ylab = "Seasonal", main = "(ii) Seasonal Component")
    abline(h = 0, col = "blue")
    
    
    plot(decomp_log_drivers$random, type = "l", col = "black", ylab = "Irregular", main = "(iii) Irregular Component")
    abline(h = 0, col = "blue") 
    
    par(mfrow = c(1, 1)) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ✅ 载入必要库
    library(dlm)
    library(dplyr)
    library(zoo)
    library(ggplot2)
    library(tidyr)
    
    # ✅ RW Kalman 模型（Local Level）+ CAPM 回归项
    estimate_rw_kalman <- function(file_path) {
      df <- read.csv(file_path) |> na.omit()
      df$Date <- as.Date(df$Date)
      y <- df$Excess_Return
      x <- df$Market_Excess_Return
      
      # 构建 build 函数
      build_rw <- function(par) {
        dlmModReg(x, dV = exp(par[2]), dW = exp(par[1]), addInt = FALSE)
      }
      
      # 固定初始值 log(1.001)
      init <- log(c(1.001, 1.001))
      fit <- dlmMLE(y, parm = init, build = build_rw, method = "BFGS")
      
      final_model <- build_rw(fit$par)
      smoothed <- dlmSmooth(y, mod = final_model)
      df$Beta_RW <- drop(smoothed$s)[-1]
      return(df)
    }
    
    # ✅ MR Kalman 模型（Mean Reverting, 正统设定）
    estimate_mr_kalman <- function(file_path) {
      df <- read.csv(file_path) |> na.omit()
      df$Date <- as.Date(df$Date)
      y <- df$Excess_Return
      x <- df$Market_Excess_Return
      n <- length(y)
      
      # 参数变换：phi 限制在 (0,1)
      transform_phi <- function(psi) (psi^2) / (1 + psi^2)
      
      build_mr <- function(par) {
        phi <- transform_phi(par[1])  # phi ∈ (0, 1)
        q1 <- exp(par[2])             # ηₜ 方差
        r <- exp(par[3])              # observation noise
        
        FF <- matrix(0, nrow = 1, ncol = 2)  # 回归项由 X 提供
        GG <- matrix(c(phi, 0, 0, 1), nrow = 2)  # 状态转移
        W <- diag(c(q1, 0))   # ζₜ ≡ 0，固定均值（MR 模型核心）
        V <- matrix(r)        # observation noise
        m0 <- c(0, 0)         # 初始状态
        C0 <- diag(1e7, 2)    # 非信息化先验
        
        mod <- dlm(FF = FF, V = V, GG = GG, W = W, m0 = m0, C0 = C0)
        mod$X <- cbind(x, rep(1, n))         # CAPM 回归项
        mod$JFF <- matrix(1, nrow = 1, ncol = 2)  # FF = β xₜ + α
        return(mod)
      }
      
      # 初始参数：φ, log(q1), log(r)
      init <- c(0, log(var(y) * 0.5), log(var(y)))
      fit <- dlmMLE(y, parm = init, build = build_mr, method = "BFGS")
      
      final_model <- build_mr(fit$par)
      smoothed <- dlmSmooth(y, final_model)
      
      df$Beta_MR <- rowSums(smoothed$s)[-1]  # 平滑估计 βₜ + α
      return(df)
    }
    
    
    # ✅ 静态 OLS beta
    theta_ols_beta <- function(industry_df, market_df) {
      df <- merge(industry_df[, c("Date", "Excess_Return")],
                  market_df[, c("Date", "Market_Excess_Return")], by = "Date")
      x <- df$Market_Excess_Return
      y <- df$Excess_Return
      sum(x * y) / sum(x^2)
    }
    
    # ✅ Rolling OLS beta
    estimate_rolling_beta <- function(industry_df, market_df, window = 52) {
      df <- merge(industry_df[, c("Date", "Excess_Return")],
                  market_df[, c("Date", "Market_Excess_Return")], by = "Date")
      df <- df[complete.cases(df), ]
      betas <- zoo::rollapply(
        1:(nrow(df) - window + 1), width = 1,
        FUN = function(i) {
          x <- df$Market_Excess_Return[i:(i + window - 1)]
          y <- df$Excess_Return[i:(i + window - 1)]
          coef(lm(y ~ x))[2]
        },
        by = 1, align = "right", fill = NA
      )
      return(data.frame(Date = df$Date, Rolling_OLS_Beta = c(rep(NA, window - 1), betas)))
    }
    
    # ✅ 绘图函数
    plot_beta_timeseries <- function(beta_df) {
      if (all(is.na(beta_df$Beta_RW))) stop("Beta_RW is all NA.")
      p <- ggplot(beta_df, aes(x = Date)) +
        geom_line(aes(y = Rolling_OLS_Beta, color = "Rolling OLS")) +
        geom_line(aes(y = Beta_MR, color = "MR Kalman")) +
        geom_line(aes(y = Beta_RW), color = "black", size = 0.8) +  # 黑色线
        geom_hline(aes(yintercept = OLS_Beta, linetype = "Static OLS"), color = "gray") +
        labs(title = "Beta Comparison", y = "Beta", color = "Model", linetype = "") +
        theme_minimal()
      if (!all(is.na(beta_df$Beta_RW))) {
        p <- p + ylim(0, 2)
      }
      print(p)
    }
    
    plot_beta_boxplot <- function(beta_df) {
      df_trimmed <- beta_df[-(1:50), ]
      box_data <- df_trimmed %>%
        select(Rolling_OLS_Beta, Beta_RW, Beta_MR) %>%
        pivot_longer(cols = everything(), names_to = "Model", values_to = "Beta")
      
      p <- ggplot(box_data, aes(x = Model, y = Beta, fill = Model)) +
        geom_boxplot() +
        geom_hline(yintercept = beta_df$OLS_Beta[1], linetype = "dashed", color = "gray") +
        labs(title = "Boxplot of Conditional Betas", y = "Beta") +
        theme_minimal()
      print(p)
    }
    
    print_beta_summary <- function(beta_df) {
      df_trimmed <- beta_df[-(1:50), ]
      summary_df <- data.frame(
        Model = c("Rolling OLS", "RW Kalman", "MR Kalman"),
        Mean = c(mean(df_trimmed$Rolling_OLS_Beta, na.rm = TRUE),
                 mean(df_trimmed$Beta_RW, na.rm = TRUE),
                 mean(df_trimmed$Beta_MR, na.rm = TRUE)),
        Range = c(diff(range(df_trimmed$Rolling_OLS_Beta, na.rm = TRUE)),
                  diff(range(df_trimmed$Beta_RW, na.rm = TRUE)),
                  diff(range(df_trimmed$Beta_MR, na.rm = TRUE))),
        Std_Dev = c(sd(df_trimmed$Rolling_OLS_Beta, na.rm = TRUE),
                    sd(df_trimmed$Beta_RW, na.rm = TRUE),
                    sd(df_trimmed$Beta_MR, na.rm = TRUE))
      )
      print(summary_df)
    }
    
    # ✅ 主流程（注意替换文件路径）
    industry_df <- read.csv("1.csv")
    industry_df$Date <- as.Date(industry_df$Date)
    market_df <- industry_df[, c("Date", "Market_Excess_Return")]
    
    ols_beta <- theta_ols_beta(industry_df, market_df)
    rolling_df <- estimate_rolling_beta(industry_df, market_df)
    rw_df <- estimate_rw_kalman("1.csv")
    mr_df <- estimate_mr_kalman("1.csv")
    
    beta_df <- mr_df[, c("Date", "Beta_MR")] %>%
      inner_join(rw_df[, c("Date", "Beta_RW")], by = "Date") %>%
      inner_join(rolling_df, by = "Date") %>%
      mutate(OLS_Beta = ols_beta)
    
    # ✅ 输出图与摘要
    plot_beta_timeseries(beta_df)
    plot_beta_boxplot(beta_df)
    print_beta_summary(beta_df)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ✅ 载入必要库
    library(dlm)
    library(dplyr)
    library(zoo)
    library(ggplot2)
    library(tidyr)
    
    # ✅ RW Kalman 模型（Local Level）
    estimate_rw_kalman <- function(file_path) {
      df <- read.csv(file_path) |> na.omit()
      df$Date <- as.Date(df$Date)
      y <- df$Excess_Return
      x <- df$Market_Excess_Return
      
      build_rw <- function(par) {
        dlmModReg(x, dV = exp(par[2]), dW = exp(par[1]), addInt = FALSE)
      }
      
      init <- log(c(1.001, 1.001))
      fit <- dlmMLE(y, parm = init, build = build_rw, method = "BFGS")
      
      final_model <- build_rw(fit$par)
      smoothed <- dlmSmooth(y, mod = final_model)
      df$Beta_RW <- drop(smoothed$s)[-1]
      return(df)
    }
    
    # ✅ MR Kalman 模型（Mean Reverting）
    estimate_mr_kalman <- function(file_path) {

      
      # 读取数据
      df <- read.csv(file_path) |> na.omit()
      df$Date <- as.Date(df$Date)
      y <- df$Excess_Return
      x <- df$Market_Excess_Return
      n <- length(y)
      
      # φ 转换函数，将 ψ 限制在 (0,1)
      transform_phi <- function(psi) (psi^2) / (1 + psi^2)
      
      # 构建 MR Kalman 模型（严格按照文献结构）
      build_mr <- function(par) {
        phi <- transform_phi(par[1])         # 均值回复系数
        q1 <- exp(par[2])                    # 状态扰动方差（只作用于 beta_t）
        r  <- exp(par[3])                    # 观测噪声方差
        
        GG <- matrix(c(phi, 0, 0, 1), nrow = 2)      # 状态转移矩阵
        W <- diag(c(q1, 0))                          # 状态噪声协方差
        V <- matrix(r)                               # 观测噪声协方差
        FF <- matrix(0, nrow = 1, ncol = 2)          # 观测矩阵
        m0 <- c(0, 0)                                 # 初始状态均值
        C0 <- diag(1e7, 2)                            # 初始状态协方差（diffuse prior）
        
        mod <- dlm(FF = FF, V = V, GG = GG, W = W, m0 = m0, C0 = C0)
        
        # ✅ 观测矩阵动态赋值：X_t = [x_t, x_t]，实现 y_t = x_t * beta_t + epsilon
        mod$X <- cbind(x, x)
        mod$JFF <- matrix(1, nrow = 1, ncol = 2)
        
        return(mod)
      }
      
      # MLE 初值设定（标准推荐）
      init <- c(0, log(var(y) * 0.5), log(var(y)))
      
      # 最大似然估计
      fit <- dlmMLE(y, parm = init, build = build_mr, method = "BFGS")
      
      # 使用估计的参数重建模型并滤波
      final_model <- build_mr(fit$par)
      smoothed <- dlmSmooth(y, final_model)
      
      # 提取 beta_t = s1 + s2 = beta_t - beta_bar + beta_bar
      df$Beta_MR <- rowSums(smoothed$s)[-1]   # 去掉初始时刻
      
      return(df)
    }
    
    estimate_mmr_kalman <- function(file_path) {

      
      # 读取数据
      df <- read.csv(file_path) |> na.omit()
      df$Date <- as.Date(df$Date)
      y <- df$Excess_Return
      x <- df$Market_Excess_Return
      n <- length(y)
      
      # φ 转换函数，确保 φ ∈ (0,1)
      transform_phi <- function(psi) (psi^2) / (1 + psi^2)
      
      # 构建 MMR Kalman 模型（完全按文献式6.21）
      build_mmr <- function(par) {
        phi <- transform_phi(par[1])         # 均值回复参数 φ
        q1 <- exp(par[2])                    # 状态1的噪声方差（η）
        q2 <- exp(par[3])                    # 状态2的噪声方差（ζ）
        r  <- exp(par[4])                    # 观测噪声方差（ε）
        
        GG <- matrix(c(phi, 0, 0, 1), nrow = 2)  # 状态转移矩阵
        W <- diag(c(q1, q2))                    # 状态扰动协方差矩阵
        V <- matrix(r)                          # 观测噪声方差
        FF <- matrix(0, nrow = 1, ncol = 2)     # 观测矩阵占位符
        m0 <- c(0, 0)                           # 初始状态均值
        C0 <- diag(1e7, 2)                      # 初始状态协方差（diffuse prior）
        
        mod <- dlm(FF = FF, V = V, GG = GG, W = W, m0 = m0, C0 = C0)
        
        # ✅ 修正关键点：观测矩阵 X_t = [x_t, x_t]
        mod$X <- cbind(x, x)
        mod$JFF <- matrix(1, nrow = 1, ncol = 2)
        
        return(mod)
      }
      
      # 初始化参数（来自文献建议）
      init <- c(0, log(var(y) * 0.5), log(var(y) * 0.01), log(var(y)))
      
      # MLE 拟合参数
      fit <- dlmMLE(y, parm = init, build = build_mmr, method = "BFGS")
      
      # 平滑估计
      final_model <- build_mmr(fit$par)
      smoothed <- dlmSmooth(y, final_model)
      
      # 状态1 + 状态2 = beta_t
      df$Beta_MMR <- rowSums(smoothed$s)[-1]
      
      return(df)
    }
    
    # ✅ 静态 OLS beta
    theta_ols_beta <- function(industry_df, market_df) {
      df <- merge(industry_df[, c("Date", "Excess_Return")],
                  market_df[, c("Date", "Market_Excess_Return")], by = "Date")
      x <- df$Market_Excess_Return
      y <- df$Excess_Return
      sum(x * y) / sum(x^2)
    }
    
    # ✅ Rolling OLS beta
    estimate_rolling_beta <- function(industry_df, market_df, window = 52) {
      df <- merge(industry_df[, c("Date", "Excess_Return")],
                  market_df[, c("Date", "Market_Excess_Return")], by = "Date")
      df <- df[complete.cases(df), ]
      betas <- zoo::rollapply(
        1:(nrow(df) - window + 1), width = 1,
        FUN = function(i) {
          x <- df$Market_Excess_Return[i:(i + window - 1)]
          y <- df$Excess_Return[i:(i + window - 1)]
          coef(lm(y ~ x))[2]
        },
        by = 1, align = "right", fill = NA
      )
      return(data.frame(Date = df$Date, Rolling_OLS_Beta = c(rep(NA, window - 1), betas)))
    }
    
    # ✅ 主流程
    industry_df <- read.csv("1.csv")
    industry_df$Date <- as.Date(industry_df$Date)
    market_df <- industry_df[, c("Date", "Market_Excess_Return")]
    
    ols_beta <- theta_ols_beta(industry_df, market_df)
    rolling_df <- estimate_rolling_beta(industry_df, market_df)
    rw_df <- estimate_rw_kalman("1.csv")
    mr_df <- estimate_mr_kalman("1.csv")
    mmr_df <- estimate_mmr_kalman("1.csv")
    
    # ✅ 合并 beta
    data_combined <- mr_df[, c("Date", "Beta_MR")] %>%
      inner_join(rw_df[, c("Date", "Beta_RW")], by = "Date") %>%
      inner_join(rolling_df, by = "Date") %>%
      inner_join(mmr_df[, c("Date", "Beta_MMR")], by = "Date") %>%
      mutate(OLS_Beta = ols_beta)
    
    # 图 1：RW Kalman vs Static OLS
    plot_rw_vs_static <- function(beta_df) {
      ggplot(beta_df, aes(x = Date)) +
        geom_line(aes(y = Beta_RW, color = "RW Kalman")) +
        geom_hline(aes(yintercept = OLS_Beta, linetype = "Static OLS"), color = "gray") +
        labs(title = "RW Kalman vs Static OLS", y = "Beta", color = "Model", linetype = "") +
        theme_minimal()
    }
    
    # 图 2：MR Kalman vs Static OLS
    plot_mr_vs_static <- function(beta_df) {
      ggplot(beta_df, aes(x = Date)) +
        geom_line(aes(y = Beta_MR, color = "MR Kalman")) +
        geom_hline(aes(yintercept = OLS_Beta, linetype = "Static OLS"), color = "gray") +
        labs(title = "MR Kalman vs Static OLS", y = "Beta", color = "Model", linetype = "") +
        theme_minimal()
    }
    
    # 图 3：MMR Kalman vs Static OLS
    plot_mmr_vs_static <- function(beta_df) {
      ggplot(beta_df, aes(x = Date)) +
        geom_line(aes(y = Beta_MMR, color = "MMR Kalman")) +
        geom_hline(aes(yintercept = OLS_Beta, linetype = "Static OLS"), color = "gray") +
        labs(title = "MMR Kalman vs Static OLS", y = "Beta", color = "Model", linetype = "") +
        theme_minimal()
    }
    
    # 图 4：Rolling OLS vs Static OLS
    plot_rolling_vs_static <- function(beta_df) {
      ggplot(beta_df, aes(x = Date)) +
        geom_line(aes(y = Rolling_OLS_Beta, color = "Rolling OLS")) +
        geom_hline(aes(yintercept = OLS_Beta, linetype = "Static OLS"), color = "gray") +
        labs(title = "Rolling vs Static OLS", y = "Beta", color = "Model", linetype = "") +
        theme_minimal()
    }
    plot_rw_vs_static(data_combined)
    plot_mr_vs_static(data_combined)
    plot_mmr_vs_static(data_combined)
    plot_rolling_vs_static(data_combined)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # 加载检验需要的包
    library(tseries)    # JB test
    library(FinTS)      # ARCH LM test
    library(lmtest)     # Ljung-Box Q test
    
    # 定义摘要函数：对 Kalman 模型提取统计量
    summarise_kalman_model <- function(model_name, fit, y, smoothed, fitted_values = NULL) {
      # 获取估计参数
      par <- fit$par
      logLik <- -fit$value
      n <- length(y)
      k <- length(par)
      bic <- log(n) * k - 2 * logLik
      
      # 提取平滑后的状态分量（假设是二维状态）
      beta <- drop(smoothed$s[, 1] + smoothed$s[, 2])
      beta_bar <- mean(beta, na.rm = TRUE)
      
      # 如果没有提供 fitted.values，就用 beta 计算拟合值
      if (is.null(fitted_values)) {
        fitted_values <- beta * y * 0 + mean(y)  # 占位；建议你传入真实的 fitted
      }
      
      # 计算 R² 拟合优度
      ss_res <- sum((y - fitted_values)^2)
      ss_tot <- sum((y - mean(y))^2)
      r2 <- 1 - ss_res / ss_tot
      
      # 参数反变换
      phi <- (par[1]^2) / (1 + par[1]^2)
      sigma_eta <- exp(par[2])
      sigma_zeta <- if (k == 4) exp(par[3]) else NA
      sigma_eps <- exp(par[k])
      
      # 各类检验
      jb <- jarque.bera.test(y)$statistic
      q12 <- Box.test(y, lag = 12, type = "Ljung-Box")$statistic
      lm6 <- ArchTest(y, lags = 6)$statistic
      
      # 汇总输出
      result <- data.frame(
        Model = model_name,
        logL = round(logLik, 2),
        BIC = round(bic, 2),
        R2 = round(r2, 3),
        phi = round(phi, 3),
        beta_bar = round(beta_bar, 3),
        `σ²_η` = round(sigma_eta, 4),
        `σ²_ζ` = if (!is.na(sigma_zeta)) round(sigma_zeta, 4) else NA,
        `σ²_ε` = round(sigma_eps, 4),
        JB = round(jb, 2),
        Q12 = round(q12, 2),
        LM6 = round(lm6, 2)
      )
      return(result)
    }
    
    
    
    
    
    
    
    
    