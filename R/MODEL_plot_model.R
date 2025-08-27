
#' Plot BMD (Benchmark Dose) Model
#'
#' This function generates a plot for a Benchmark Dose (BMD) model.
#'
#' @param model The fitted BMD model.
#' @param cex The size of points in the plot (default: 6).
#' @param x_pos The position of the x-axis (default: 15).
#' @param y_pos_th The position of the y-axis threshold (default: 0.85).
#' @param confidence_interval The confidence interval level (default: 0.95).
#' @param title_label The title of the plot (default: "title").
#' @return A plot displaying the BMD model with confidence intervals and other related data points.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
plot_bmdx = function(model,
         cex = 6,
         confidence_interval = 0.95,
         plot_ic50 = FALSE,
         xlim_u = NULL,
         title_label = "title"){

  data = data.frame(dose = model$data_x, expr = model$data_y)
  x_breaks = sort(unique(data$dose))
  new.data <- data.frame(dose = seq(min(model$data_x), max(model$data_x), length.out = 1000))

  interval = as_tibble(predict(model, newdata = new.data, interval = "confidence",level = confidence_interval)) %>%
    mutate(dose = new.data$dose)

  p1 <- ggplot2::ggplot(data) +
    ggplot2::geom_point(ggplot2::aes(x = dose, y = expr),size = 2) +
    ggplot2::xlab("Concentration/Dose") +
    ggplot2::ylab("Expr")

  if(!is.null(xlim_u)){
    p1 = p1+xlim(0, xlim_u)
  }

  p1 = p1 + ggplot2::geom_line(data     = interval, ggplot2::aes(x = dose, y = fit )) +
    ggplot2::geom_ribbon(data   = interval, ggplot2::aes(x = dose, ymin = lwr, ymax = upr), alpha = 0.5, inherit.aes = F, fill = "gray")

  if (is.null(model$BMD) == FALSE) {
    geneBMD = model$BMD
    geneBMDL = model$BMDL
    geneBMDU = model$BMDU
    icval = model$AC50

    x_breaks = c(x_breaks,round(c(geneBMD, geneBMDL, geneBMDU, icval),2))

    pred_bmd  = predict(model, newdata = data.frame(dose = geneBMD))
    pred_bmdu = predict(model, newdata = data.frame(dose = geneBMDU))
    pred_ic50 = predict(model, newdata = data.frame(dose = icval))

    min_expr_value = min(data$expr)

    if(plot_ic50){
      bmd_df = data.frame(x     = c(geneBMD, geneBMDL, geneBMDU, 0, icval, 0),
                          y     = c(rep(min_expr_value,3),pred_bmd, min_expr_value,pred_ic50),
                          xend  = c(geneBMD, geneBMDL, geneBMDU,geneBMDU, icval, icval),
                          yend  = c(rep(pred_bmd,4), pred_ic50, pred_ic50),
                          label = c("BMD", "BMDL", "BMDU", "BMR", "IC50/EC50","IC50/EC50"))
    }else{
      bmd_df = data.frame(x     = c(geneBMD, geneBMDL, geneBMDU, 0),
                          y     = c(rep(min_expr_value,3),pred_bmd),
                          xend  = c(geneBMD, geneBMDL, geneBMDU,geneBMDU),
                          yend  = c(rep(pred_bmd,4)),
                          label = c("BMD", "BMDL", "BMDU", "BMR"))
    }

    if(is.na(geneBMD) & is.na(geneBMDL) & is.na(geneBMD)){
      p1 = p1
    }else{
      if(plot_ic50){
        p1 =  p1 +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data = bmd_df, linetype = "dashed", lwd = 1) +
          geom_point(aes(x = geneBMD, y = pred_bmd), size = 2, color = "blue") +
          geom_point(aes(x = geneBMDL, y = pred_bmd), size = 2, color = "red") +
          geom_point(aes(x = geneBMDU, y = pred_bmd), size = 2, color = "green") +
          geom_point(aes(x = icval, y = pred_ic50), size = 2, color = "orange")
      }else{
        p1 =  p1 +
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data = bmd_df, linetype = "dashed", lwd = 1) +
          geom_point(aes(x = geneBMD, y = pred_bmd, color = "blue"), size = 2, color = "blue") +
          geom_point(aes(x = geneBMDL, y = pred_bmd, color = "red"), size = 2, color = "red") +
          geom_point(aes(x = geneBMDU, y = pred_bmd, color ="green"), size = 2, color = "green")
      }

    }


    if(plot_ic50){
      # p1 = p1 + scale_color_manual(values = c(Blue = "blue", Red = "red", Green = "green", Orange = "orange"),
      #                              labels = c("BMD", "BMDL", "BMDU", "IC50"))

      p1 = p1 + scale_color_manual(name = "Legend",
                                   values = c("BMD" = "blue",
                                              "BMDL" = "red",
                                              "BMDU" = "green",
                                              "IC50" = "orange"))
    }else{
      # p1 = p1 + scale_color_manual(values = c(Blue = "blue", Red = "red", Green = "green"),
      #                              labels = c("BMD", "BMDL", "BMDU"))

      p1 = p1 + scale_color_manual(name = "Legend",
                         values = c("BMD" = "blue",
                                    "BMDL" = "red",
                                    "BMDU" = "green"))

    }

    p1 = p1 + ggplot2::theme(plot.title   = ggplot2::element_text(size = cex, face = "bold"),
                   legend.title = ggplot2::element_text(size = cex),
                   legend.text  = ggplot2::element_text(size = cex),
                   axis.text    = ggplot2::element_text(size = cex),
                   axis.title   = ggplot2::element_text(size = cex,face = "bold"),
                   legend.position = "right") +
      ggtitle(title_label)




    label = paste0("BMD = ", round(geneBMD,2), "\nBMDL = ", round(geneBMDL,2),"\nBMDU = ", round(geneBMDU,2), "\nIC50/EC50 = " , round(icval,2), sep = "")
  }
  p1
}

# plot_bmdx = function(model,
#                      cex = 6,
#                      x_pos = 15,
#                      y_pos_th = 0.85,
#                      confidence_interval = 0.95,
#                      title_label = "title"){
#
#   data = data.frame(dose = model$data_x, expr = model$data_y)
#   x_breaks = sort(unique(data$dose))
#   new.data <- data.frame(dose = seq(min(model$data_x), max(model$data_x), length.out = 1000))
#
#   interval = as_tibble(predict(model, newdata = new.data, interval = "confidence",level = confidence_interval)) %>%
#     mutate(dose = new.data$dose)
#
#   p1 <- ggplot2::ggplot(data) +
#         ggplot2::geom_point(ggplot2::aes(x = dose, y = expr),size = 2, colour = "black") +
#         ggplot2::xlab("Concentration/Dose") +
#         ggplot2::ylab("Expr")
#
#   p1 = p1 + ggplot2::geom_line(data     = interval, ggplot2::aes(x = dose, y = fit )) +
#             ggplot2::geom_ribbon(data   = interval, ggplot2::aes(x = dose, ymin = lwr, ymax = upr), alpha = 0.5, inherit.aes = F, fill = "gray") +
#             ggplot2::theme(plot.title   = ggplot2::element_text(size = cex, face = "bold"),
#                            legend.title = ggplot2::element_text(size = cex),
#                            legend.text  = ggplot2::element_text(size = cex),
#                            axis.text    = ggplot2::element_text(size = cex),
#                            axis.title   = ggplot2::element_text(size = cex,face = "bold"), legend.position = "none")
#
#   if (is.null(model$BMD) == FALSE) {
#     geneBMD = model$BMD
#     geneBMDL = model$BMDL
#     geneBMDU = model$BMDU
#     icval = model$AC50
#
#     x_breaks = c(x_breaks,round(c(geneBMD, geneBMDL, geneBMDU, icval),2))
#
#     pred_bmd  = predict(model, newdata = data.frame(dose = geneBMD))
#     pred_bmdu = predict(model, newdata = data.frame(dose = geneBMDU))
#     pred_ic50 = predict(model, newdata = data.frame(dose = icval))
#
#     min_expr_value = min(data$expr)
#
#     bmd_df = data.frame(x     = c(geneBMD, geneBMDL, geneBMDU, 0, icval, 0),
#                         y     = c(rep(min_expr_value,3),pred_bmd, min_expr_value,pred_ic50),
#                         xend  = c(geneBMD, geneBMDL, geneBMDU,geneBMDU, icval, icval),
#                         yend  = c(rep(pred_bmd,4), pred_ic50, pred_ic50),
#                         label = c("BMD", "BMDL", "BMDU", "BMR", "IC50/EC50","IC50/EC50"))
#
#     p1 = p1 + ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend),data = bmd_df,linetype = "dashed", lwd = 1) +
#               ggplot2::geom_point(  ggplot2::aes(x = geneBMD, y = pred_bmd), size = 2,col = "black") +
#               ggplot2::geom_point(  ggplot2::aes(x = geneBMDL, y = pred_bmd), size = 2,col = "black") +
#               ggplot2::geom_point(  ggplot2::aes(x = geneBMDU, y = pred_bmd), size = 2,col = "black") +
#               ggplot2::scale_color_manual(name = "",values = c("red","blue","green","black"),
#                                           breaks = c( "BMDL","BMD","BMDU","BMR"),
#                                           labels = c( "BMDL","BMD","BMDU","BMR")) + ggtitle(title_label)
#
#     label = paste0("BMD = ", round(geneBMD,2), "\nBMDL = ", round(geneBMDL,2),"\nBMDU = ", round(geneBMDU,2), "\nIC50/EC50 = " , round(icval,2), sep = "")
#   }
#   p1
# }
