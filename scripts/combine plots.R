## utility script to combine various figures

library(magick)


julian.2sg <- image_read("figs/julian_predicted_effect_cod2sg_zinb.png")
ssb.2sg <- image_read("figs/SSB_predicted_effect_cod2sg_zinb.png")
temp.2sg <- image_read("figs/temp.anom_predicted_effect_cod2sg_zinb.png")
img <- c(julian.2sg, ssb.2sg, temp.2sg)

stack1 <- image_append(image_scale(img, "100%"))
stack1 <- image_annotate(stack1, "k=5", color = "red", size = 70)

julian.2sg3 <- image_read("figs/julian_predicted_effect_cod2sg3_zinb.png")
ssb.2sg3 <- image_read("figs/SSB_predicted_effect_cod2sg3_zinb.png")
temp.2sg3 <- image_read("figs/temp.anom_predicted_effect_cod2sg3_zinb.png")
img <- c(julian.2sg3, ssb.2sg3, temp.2sg3)

stack2 <- image_append(image_scale(img, "100%"))
stack2 <- image_annotate(stack2, "k=3", color = "red", size = 70)
img <- c(stack1, stack2)
stack <- image_append(image_scale(img, "100%"), stack = T)

image_write(stack, path = "figs/cod2sg_cod2sg3_plot_comparison.png", format = "png")


## combine temp.anom and cod model recruitment vs. FAR plots ------------------------

temp.anom <- image_read("./figs/temp.anom_predicted_effect_cod2sg_zinb_k3.png")
FAR <- image_read("./figs/prelim_FAR_recruit_plot.png")

img <- c(temp.anom, FAR)

stack1 <- image_append(image_scale(img, "100%"))

R.FAR <- image_read("./figs/predicted_effect_cod_R_FAR.png")
stack2 <- image_append(image_scale(R.FAR, "100%"))
img <- c(stack1, stack2)
stack <- image_append(image_scale(img, "100%"), stack = T)
image_write(stack, path = "figs/cod_R_FAR_stack.png", format = "png")


## combine cod-pollock R plot, pollock seine vs. FAR and pollock model recruitment vs. FAR plots ------------

R.R <- image_read("./figs/predicted_effect_cod_poll_R.png")
FAR <- image_read("./figs/pollock_FAR_recruit_plot.png")

img <- c(R.R, FAR)

stack1 <- image_append(image_scale(img, "100%"))

R.FAR <- image_read("./figs/predicted_effect_pollock_R_FAR.png")
stack2 <- image_append(image_scale(R.FAR, "100%"))
img <- c(stack1, stack2)
stack <- image_append(image_scale(img, "100%"), stack = T)
image_write(stack, path = "figs/pollock_R_FAR_stack.png", format = "png")


## combine modeled CMIP FAR projections and cod-pollock projected R plot ------------------------
plot.nil <- ggplot() + theme_void()

png("./figs/Fig4-projected_FARandR.png", width=5, height=5, units='in', res=300)
ggpubr::ggarrange(CMIP.FAR, 
                  ggpubr::ggarrange(cod.poll.proj.R, plot.nil, widths=c(0.7, 0.3)),
                            ncol=1, heights=c(1,0.9))
dev.off()
