
  boxPlotSingle =function(data, main = ''){
    #Aggregate by ccs level
    means = aggregate(score_bm~CCSLEVEL, data,mean)

    means[ ,2] = round(means[ ,2],2)
    p = ggplot(data, aes(x=CCSLEVEL, y = score_bm, color = CCSLEVEL)) +
       geom_boxplot() +   ggtitle(main) +
     geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, binwidth = .25) +
     stat_summary(fun.y=mean, colour="darkred", geom="point",
                    shape=18, size=3) +
       geom_text(data = means, aes(label = score_bm, y = score_bm + 0.15))

    p
  }




visualizeReturner = function(returner, mat = 'inst_mat',
                  type = 'score', only_significant = FALSE,
                  plot_order = 'effect', xlab = '', ylab = '',lwd = 1,main = '', ...){ #given/ranked, some other other order, 'sorted'

 #Extract id colums for each option
  if (mat == 'inst_mat'){
    id_cols = 'inst'
  }	 else if (mat == 'subset_mat_baseline'){
    id_cols = 'subset'
  } else if (mat == 'subset_mat_nobaseline'){
    id_cols = 'subset'
  } else if ( mat == 'inst_subset_mat_baseline'){
    id_cols = c('inst', 'subset')
  } else if (mat == 'inst_subset_mat_nobaseline'){
    id_cols = c('inst', 'subset')
  }
  if (type == 'score'){
    cols = c('score_est', 'score_lower', 'score_upper')
  } else if (type == 'effect'){
    cols = c('effect_est', 'effect_lower', 'effect_upper')
  }

  mat_use = returner[[mat]]
  id_vec = mat_use[ ,id_cols]
  n = dim(mat_use)[1]
  if (length(id_cols) == 2){
    id_vec = apply(id_vec, 1, idStr)
  }

  #ID vars
  plot_vars = mat_use[ ,cols] #Subse

  #First order, then subset
  if (plot_order == 'ranked' ){
    order_vec = order(plot_vars[ ,1])
  } else if (plot_order == 'source'){
	order_vec = 1:n
  } else if (plot_order == 'effect'){
	order_vec = order(mat_use$effect_est)
  }

  indices = is.finite(plot_vars[ ,1])
  if (only_significant){
	indices = plot_vars[ ,2] > 0 | plot_vars[ ,3] < 0
  }
    #Extract the variables we want at the end
	plot_vars = plot_vars[order_vec, ]
	plot_names = id_vec[order_vec]
	indices = indices[order_vec]
	plot_vars = plot_vars[indices, ]
	plot_names = plot_names[indices]
  center = plot_vars[ ,1]
  lower = plot_vars[ ,2]
  upper = plot_vars[ ,3]
  p = length(center)
  x = 1:p

  graphics::plot(center, ylim = c(min(lower), max(upper)),
  xlab = xlab, ylab = '', axes = FALSE)
  box()
  axis(2)
  axis(1, at = x, labels=plot_names, las=2, cex.axis = 0.65)
  mtext(side = 2, line = 2, ylab)
  graphics::segments(x0=x, x1=x,y0 = lower, y1 = upper, lwd = lwd, ...)
  lower_indices = upper < 0
  upper_indices = lower > 0
  graphics::segments(x0=x[lower_indices], x1=x[lower_indices],
                     y0 = lower[lower_indices], y1 = upper[lower_indices],  col = 'red',lwd = lwd* 1.5 , ...)
  graphics::segments(x0=x[upper_indices], x1=x[upper_indices],
                     y0 = lower[upper_indices], y1 = upper[upper_indices],  col = 'blue', lwd = lwd * 1.5, ...)
  graphics::abline(h = 0, lty = 3)
  graphics::title(main)

  #Order and subset? if necesasary
  #Here is where we would extrreact only_significatnt
  return(0)
}

# comparePlot = function(x, y1, y2=NULL, col_vec = c('blue', 'grey'),
                       # rank = TRUE,
                       # xlab = '', ylab1 = '', ylab2 = '', main = '', pch = 19, plot_legend = TRUE, ...){
  # if (rank){
			# x=rank(x); y1=rank(y1); y2=rank(y2)
			# }
  # plot(x, y1, xlab = '', ylab = '', main = main , pch = pch, col = col_vec[1])
  # mtext(ylab1, side = 2, line = 1.6)
  # mtext(ylab2, side = 4, line = 0)
  # legend_add = NULL
  # if (length(y2) > 0){points(x, y2, col = col_vec[2], pch = pch)
  # legend_add = paste0(ylab2, ': r=', round(cor(x,y2),3)) }
  # points(x, y1, col = col_vec[1], pch = pch)
  
  # if (plot_legend){
  # #Build legend manually
  # legend_str = c(
	# paste0(ylab1, ': r=', round(cor(x,y1),3)),
	# legend_add)
  # legend('topleft', legend_str, fill = col_vec)
  # }
# }

comparePlot = function(x, y1, y2=NULL, col_vec = c('blue', 'grey'),
                       rank = TRUE,
                       xlab = '', ylab1 = '', ylab2 = '', main = '', pch = 19, plot_legend = TRUE, ...){
  if (rank){
			x=rank(x); y1=rank(y1); y2=rank(y2)
			}
			
  plot(x, y1, xlab = '', ylab = '', main = main , pch = pch, col = col_vec[1])
  mtext(ylab1, side = 2, line = 1.6)
  mtext(ylab2, side = 4, line = 0)
  legend_add = NULL
  if (length(y2) > 0){points(x, y2, col = col_vec[2], pch = pch)
  legend_add = paste0(ylab2, ': r=', round(cor(x,y2),3)) }
  points(x, y1, col = col_vec[1], pch = pch)
  
  if (plot_legend){
  #Build legend manually
  legend_str = c(
	paste0(ylab1, ': r=', round(cor(x,y1),3)),
	legend_add)
  legend('topleft', legend_str, fill = col_vec)
  }
}

comparePlotgg = function(data_frame, x = 'score_bm', y = 'score_lr', 
	main = '',
	col = 'CCSLEVEL', col_as_factor = TRUE,
	size_var = 'n'){

#Convert the color variable to a factor...
if (col_as_factor){data_frame[[col]] = factor(data_frame[[col]])}
ggplot(data_frame, aes_string(x = x, y = y)) + geom_point(aes_string(color = col)) + ggtitle(main)

}