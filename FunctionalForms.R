################################################################################
#                      Computational Model of Acrophobia                       #
#                 Visualizing the relations between variables                  #
################################################################################ 


#--- Set up ----------------

text_size = 2
title_size = 3

pdf("Fig3.pdf", height = 11, width = 12)
par(mfrow = c(2,2))
par(mar = c(8, 7, 4, 2)) # increase left margin 

#--------------------------

# =========================
# Meters on Punishment value 
# =========================

punishment_h <- 0.5
meters <- seq(0, 100, 1)
punishment <- -1 / (1 + (exp(1)^(-5 * (log10(meters) - 2 * punishment_h))))

plot(meters, punishment, type = "l", xlab = "Height (in meters)", ylab = "Punishment value", axes = FALSE, cex.lab = text_size, ylim = c(-1,0))
axis(side = 1, cex.axis = text_size); axis(side = 2, cex.axis = text_size)

# add different parameter value:
punishment_h <- 0.3
punishment <- -1 / (1 + (exp(1)^(-5 * (log10(meters) - 2 * punishment_h))))
lines(meters, punishment, type = "l", lty = 2) 

# add different parameter value: 
punishment_h <- 0.7
punishment <- -1 / (1 + (exp(1)^(-5 * (log10(meters) - 2 * punishment_h))))
lines(meters, punishment, type = "l", lty = 3) 

# add legend: 
legend('topright', legend = c('h = 0.3', 'h = 0.5', 'h = 0.7'),
       lty = c(2, 1, 3), bty = 'n', cex = text_size)
title(main = "(A)", cex.main = title_size)

# =========================
# Q^p(s,a) on Perceived Threat 
# =========================

Q <- seq(-1, 0, .001)
plot(Q, abs(Q), type = "l", xlab = "Negative Q-value", ylab = "Perceived Threat", axes = FALSE, cex.lab = text_size)
axis(side = 1, cex.axis = text_size); axis(side = 2, cex.axis = text_size)
title(main = "(B)", cex.main = title_size)

# =========================
# Perceived Threat on Bodily Responses 
# =========================

T <- seq(0, 1, .001) 
plot(T, T^.1, type = "l", xlab = "Perceived Threat", ylab = "Change in Bodily Responses", axes = FALSE, cex.lab = text_size)
axis(side = 1, cex.axis = text_size); axis(side = 2, cex.axis = text_size)
lines(T, T^.25, type = "l", lty = 2) 
lines(T, T^.4, type = "l", lty = 3)
legend('bottomright', legend = c('a = 0.1', 'a = 0.25', 'a = 0.4'),
       lty = c(1, 2, 3), bty = 'n', cex = text_size)
title(main = "(C)", cex.main = title_size)

# =========================
# Arousal on Sensation
# =========================

# A <- seq(0, 1, 0.001) 
# S <- 1 / (1 + exp(-10*(A-.5)))
# 
# plot(A, S, type = "l", xlab = "Bodily Sensations", ylab = "Change in Somatic Sensation", axes = FALSE, cex.lab = text_size)
# axis(side = 1, cex.axis = text_size); axis(side = 2, cex.axis = text_size)
# title(main = "(D)", cex.main = title_size)

# =========================
# Bodily Responses on Perceived Control
# =========================

C_max <- .9
S <- seq(0, 1, 0.001)
C <- C_max / (1 + exp(10 * (S - 0.5)))

plot(S, C, type = "l", xlab = "Bodily Responses", ylab = "Change in Perceived Control", axes = FALSE, cex.lab = text_size, ylim = c(0,1))
axis(side = 1, cex.axis = text_size); axis(side = 2, cex.axis = text_size)
lines(S, rep(C_max, length(S)), lty = 2)
text(0.9, C_max + 0.02, "C max", pos = 3, cex = text_size) # add text above C_max line
title(main = "(E)", cex.main = title_size)


# =========================
# Fear on Perceived Reward 
# =========================

# F <- seq(0, 1, 0.001)
# 
# plot(F, F^0.5, type = "l", xlab = "Fear", ylab = "Perceived Reward", axes = FALSE, cex.lab = text_size)
# axis(side = 1, cex.axis = text_size); axis(side = 2, cex.axis = text_size)
# lines(F, F, type = "l", lty = 2) 
# lines(F, F^2, type = "l", lty = 3)
# legend('bottomright', legend = c('\u03B1 < 1', '\u03B1 = 1', '\u03B1 > 1'),
#       lty = c(1, 2, 3), bty = 'n', cex = text_size)
# title(main = "(F)", cex.main = title_size)

#--------------------------

dev.off()

