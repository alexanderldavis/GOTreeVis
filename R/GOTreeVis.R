#' Tree visualization for Gene Ontology enrichment results
#'
#' \code{GOTreeVis} creates a GOTree visualisation of GO enrichment data and writes it to a pdf file.
#'
#' @param disp_data Enrichment data to visualize.
#' Four columns: \enumerate{
#' \item Term
#' \item Pvalue - scales the length of the branches
#' \item Count (i.e. number of genes matching the term) - scales the thickness of the branches.
#' If this scaling is not needed, this column should be set to all equal values. All branches will then be set to min_size.
#' \item a column used to determine the side of the tree, e.g. ontology (BP/MF).
#' }
#' In case of a one-sided tree (see oneside parameter) the fourth column can be omitted.
#' The side column should be a factor with two levels (first level = left, second = right).
#' For character side columns a factor will be generated (levels will be in alphabetical order).
#' @param out_file Filename of the output file. The name must end with a valid file extension to define the graphics device. PDF output is recommended.
#' @param b_rootlabel Text for the label at the root of the tree. [default: ""]
#' @param b_sidelabels Text for the labels on each side; The first value corresponds to the left side. [default: ('GO:MF','GO:BP')]
#' @param b_parse Parse text for root and side labels using plotmath. Can be used to plot labels containing superscript or subscript.
#' @param root_color    Color for the 'root' box at the bottom of the tree. Applies to text and border. [default: black]
#' @param colors        Colors for the two sides of the tree, again first=left. [default: blue/green]
#' @param plot_pvalue_labels  Should the p-value labels on the axis be plotted or not? [default: T]
#' @param font          Font family to be used for the text in the plot [default: sans]
#' @param oneside       Special functionality in case only one category should be plotted.
# possible values: 'l' or 'r', meaning all data will be plotted left or right.
#' @param pval_scaling      Scaling factor between -log(pvalue) and branch length; Gets generated automatically from the dimensions of the data if left unset (recommended for first run). Can be used to scale the tree in x dimension. [default: NA]
#' @param br_spacing        Spacing between branches. Gets generated automatically from the dimensions of the data if left unset (recommended for first run). Can be used to scale the tree in y dimension. [default: NA]
#' @param tr_width          Adjusts the width of the trunk. [default: 4]
#' @param br_min_size       Sets the count->branch_width scale. Defines the size(width) of the branch with the smallest count value. [default: 2]
#' @param br_max_size       Sets the count->branch_width scale. Defines the size(width) of the branch with the highest count value. [default: 4]
#' @param br_angle          Angle of the branches in radians as multiples of pi. [default: 1/7, equal to 180°/7]
#' @param br_textoffset_x   Manually offset the text from the branches (x-axis). [default: 0]
#' @param br_textoffset_y   Manually offset the text from the branches (y-axis). [default: 0]
#' @param br_offset         Adjust the position of the top branch (e.g. if the position should be at the very top). [default: scales with spacing]
#' @param plot_size         Plot size, passed to ggsave for both height and width; changes of this value will likely require adjustments of label_size. [default: 8]
#' @param label_size        Adjusts the font size of all labels in the plot. [default: 4]
#' @param expand_x          [default: (0.1,0.1)] can be used to expand the plot area on each side in case of long labels
#'
#' @return PDF file (out_file)
#' @import ggplot2
#' @export
GOTreeVis <- function(disp_data, out_file,
                         b_rootlabel = "", b_sidelabels = c("GO:MF", "GO:BP"), b_parse = F,
                         plot_pvalue_labels = T,root_color = "black", colors = c("deepskyblue3", "seagreen3"),
                         oneside = NULL, font = "sans",
                         #pval_scaling = 1.9, br_spacing = 5,
                         pval_scaling = NA, br_spacing = NA,
                         tr_width = 4, br_min_size = 2, br_max_size = 4,
                         br_angle = 1/7, br_textoffset_x = 0, br_textoffset_y = 0, br_offset = NA,
                         plot_size = 8, label_size = 4, expand_x = c(0.1, 0.1)) {
    # internal parameters
    tr_x <- 0
    tr_yst <- 0
    text_col <- 1
    pval_col <- 2
    count_col <- 3
    side_col <- 4

    # set up side column, especially in case of a onesided plot.  in case of onesided plot expand_x is changed (plot is less wide, so
    # higher numbers are needed to expand)
    if (!is.null(oneside)) {
      disp_data[[side_col]] <- factor(rep(oneside, nrow(disp_data)), levels = c("l", "r"))
      if(oneside == "r") expand_x <- c(0.2, 0.6)
      else expand_x <- c(0.6, 0.2)
    }

    # check correctness of data.frame to display.
    if (ncol(disp_data) != 4 |
        ! is.character(disp_data[[text_col]]) | ! is.numeric(disp_data[[pval_col]]) |
        ! is.numeric(disp_data[[count_col]])) {
      stop("Input data.frame columns have to be: Term (character), Pvalue (numeric), Count (numeric), side_col!")
    }
    if (br_min_size >= br_max_size) {
      stop("br_min_size must be different from and smaller than br_max_size!")
    }
    if(!is.factor(disp_data[[side_col]]))
      disp_data[[side_col]] <- as.factor(disp_data[[side_col]])

    # set ranks to set branch locations
    disp_data$rank <- 0
    oldO <- disp_data[1, side_col]
    oldR <- 0
    for (i in 2:nrow(disp_data)) {
        currO <- disp_data[i, side_col]
        if (currO == oldO) {
            disp_data[i, "rank"] <- oldR + 2
        } else {
            disp_data[i, "rank"] <- oldR + 1
        }
        oldR <- disp_data[i, "rank"]
        oldO <- currO
    }
    # derive base values for spacing and p-value scaling from data
    if (is.na(br_spacing)){
      br_spacing <- round(54 / nrow(disp_data),0)
      print(paste("br_spacing set to",br_spacing))
    }

    # set offset for first branch (if not set it's slightly offset from the top)
    if (is.na(br_offset))
        br_offset <- br_spacing / 2.7

    # trunkdata - plotdata for trunk #### set trunk length based on br_spacing
    tr_len <- br_spacing * max(disp_data$rank) + br_spacing / 3
    tr_x_offs <- tr_width / 4
    tr_x_l <- tr_x - tr_x_offs
    tr_x_r <- tr_x + tr_x_offs
    trunkdata <- data.frame(xmin = c(tr_x_l, tr_x_r) - tr_x_offs, xmax = c(tr_x_l, tr_x_r) + tr_x_offs,
                            ymin = tr_yst, ymax = tr_yst + tr_len, side = c("l", "r"), stringsAsFactors = F)
    if (!is.null(oneside))
        trunkdata <- trunkdata[trunkdata$side == oneside, ]

    # base - plotdata for root label ####
    base <- data.frame(xst = tr_x, yst = tr_yst, text = b_rootlabel)

    # treedata - plotdata for branches ####

    # offset for first branch to not start at the very top
    top_branch <- tr_yst + tr_len - br_offset

    # define scale for count (translates to branch width)
    min_count <- min(disp_data[[count_col]])
    max_count <- max(disp_data[[count_col]])
    # if all values are the same, set dist to 1, in this case all will be scaled to min_size
    dist_count <- ifelse(max_count == min_count, 1, max_count - min_count)
    dist_size <- br_max_size - br_min_size

    # for text labels: offset from the trunk
    x_off <- 2.5 + tr_width/8

    # generate branch data
    treedata <- plyr::adply(disp_data, .margins = 1, .id = NULL, .expand = F, function(row) {
        yst <- top_branch - row[["rank"]] * br_spacing
        lp <- log(row[[pval_col]], base = 10)  # log10(pvalue)

        # rescale count to size scale (defined by br_min_size/br_max_size)
        size <- br_min_size + (row[[count_col]] - min_count)/dist_count * dist_size
        # scale y-offset of text_labels with branch size
        text_y_offset <- 0.7 + br_textoffset_y + size/4 + ifelse(is.null(oneside),0,0.2)

        if (row[[side_col]] == levels(row[[side_col]])[1]) {
            # left branches
            br_offset <- br_offset + br_spacing/2
            textangle = -br_angle
            lineangle = 1 - br_angle
            just = 1
            side = "l"
            trunk_x = tr_x_l
            text_x_offset = -(br_textoffset_x + x_off)
        } else {
            # right branches
            textangle = br_angle
            lineangle = br_angle
            just = 0
            side = "r"
            trunk_x = tr_x_r
            text_x_offset = br_textoffset_x + x_off
        }
        newBranch <- data.frame(xst = trunk_x, yst = yst,
                                angle = lineangle,textangle = textangle,
                                len = -lp, size = size,
                                text = row[[text_col]], text_x_offset = text_x_offset,
                                text_y_offset = text_y_offset, just = just,
                                side = side, stringsAsFactors = F)
        return(newBranch)
    })

    # get maximum -log10(pvalue) for each side
    if (is.null(oneside)) {
      max_pval_l <- max(treedata[treedata$side == "l", ]$len)
      max_pval_r <- max(treedata[treedata$side == "r", ]$len)
    } else {
      max_pval_l <- max(treedata$len)
      max_pval_r <- max_pval_l
    }
    # maximum p value displayed on the axis (next multiple of 5 after the maximum in the data; at least 10^-10)
    ax_pval_max_l <- max(10, ceiling(max_pval_l/5) * 5)
    ax_pval_max_r <- max(10, ceiling(max_pval_r/5) * 5)

    if(is.na(pval_scaling)) {
      sum_x <- ax_pval_max_l + ifelse(is.null(oneside), ax_pval_max_r, 0)
      pval_scaling <- round(tr_len / sum_x * 1.5,1)
      print(paste("pval_scaling set to",pval_scaling))
    }
    # scale p-value:
    treedata$len <- treedata$len * pval_scaling


    # baseaxis - plotdata for p value axis####
    # calculate x-axis projections for largest p-values to derive p-value scale on the x-axis (fact_pval_to_xlen).
    if (!is.null(oneside)) {
        if (oneside == "r") {
            graph_xst <- min(treedata$xst - treedata$len * cospi(treedata$angle))
            graph_xen <- max(treedata$xst + treedata$len * cospi(treedata$angle))
            len_l <- tr_x_l - graph_xst
            len_r <- graph_xen - (tr_x_r)
            fact_pval_to_xlen <- len_r/max_pval_r
        } else if (oneside == "l") {
            graph_xst <- min(treedata$xst + treedata$len * cospi(treedata$angle))
            graph_xen <- max(treedata$xst - treedata$len * cospi(treedata$angle))
            len_l <- tr_x_l - graph_xst
            len_r <- graph_xen - (tr_x_r)
            fact_pval_to_xlen <- len_l/max_pval_l
        }
    } else {
        graph_xst <- min(treedata$xst + treedata$len * cospi(treedata$angle))
        graph_xen <- max(treedata$xst + treedata$len * cospi(treedata$angle))
        len_l <- tr_x_l - graph_xst
        len_r <- graph_xen - (tr_x_r)
        fact_pval_to_xlen <- len_r/max_pval_r  # equivalent to len_l/max_pval_l (same scale left and right)
    }
    # define length of axis by scaling the maximum values displayed
    ax_len_l <- fact_pval_to_xlen * ax_pval_max_l
    ax_len_r <- fact_pval_to_xlen * ax_pval_max_r
    tick_x <- c(-seq(5, ax_pval_max_l, 5), seq(5, ax_pval_max_r, 5))  # create ticks every 5 orders of magnitude

    baseaxis <- data.frame(xst = c(tr_x_l - ax_len_l, tr_x_r), yst = tr_yst, angle = 0,
                           len = c(ax_len_l, ax_len_r), text = "", stringsAsFactors = F)
    # add extra line between left and right to make one continuous axis
    baseaxis <- rbind(baseaxis, data.frame(xst = tr_x_l, yst = tr_yst, angle = 0,
                                           len = tr_x_r - tr_x_l, text = "", stringsAsFactors = F))
    tick_len <- 1  # length of each tick
    # add plotdata for each tick
    for (tick in tick_x) {
        baseaxis <- rbind(baseaxis, data.frame(xst = ifelse(sign(tick) < 0, tr_x_l, tr_x_r) + tick * fact_pval_to_xlen,
                                               yst = tr_yst,
            angle = 3/2, len = tick_len, text = paste0("10^{", -abs(tick), "}"), stringsAsFactors = F))
    }

    # baselabels - plotdata for side labels####

    # position labels between ticks in the middle of the axis (same position on both axes, determined by shorter axis)
    ax_l<-baseaxis[baseaxis$text!="" & baseaxis$xst < tr_x,"xst"]
    ax_r<-baseaxis[baseaxis$text!="" & baseaxis$xst > tr_x,"xst"]
    if(length(ax_l) > length(ax_r))
      mid<-as.integer(length(ax_r)/2)+1
    else
      mid<-as.integer(length(ax_l)/2)+1
    # position between two ticks:
    labeltick_l <- ax_l[mid] - (ax_l[mid] - ax_l[mid-1])/2
    labeltick_r <- ax_r[mid] - (ax_r[mid] - ax_r[mid-1])/2
    bl_x <- c(labeltick_l, labeltick_r)
    baselabels <- data.frame(x = bl_x, y = rep(tr_yst + 0.45, 2), just = 0.5,
                             text = b_sidelabels, side = c("l", "r"), stringsAsFactors = F)
    if (!is.null(oneside)) {
        if (oneside == "l") {
            baseaxis <- baseaxis[baseaxis$xst < tr_x_l, ]
            baselabels <- baselabels[baselabels$x < tr_x_l, ]
        } else if (oneside == "r") {
            baseaxis <- baseaxis[baseaxis$xst > tr_x_l, ]
            baselabels <- baselabels[baselabels$x > tr_x_l, ]
        }
    }

    # create factors here to ensure they have both levels and colors get assigned correctly (especially relevant for the oneside case)
    trunkdata$side <- factor(trunkdata$side, levels = c("l", "r"))
    treedata$side <- factor(treedata$side, levels = c("l", "r"))
    baselabels$side <- factor(baselabels$side, levels = c("l", "r"))

    # plot data #### draw the tree trunk
    gg <- ggplot() + geom_rect(data = trunkdata, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = side))
    # draw the tree branches
    gg <- gg + geom_segment(data = treedata[nrow(treedata):1, ], aes(x = xst, y = yst,
                                                                     xend = xst + cospi(angle) * len,
                                                                     yend = yst + sinpi(angle) * len,
                                                                     size = size, color = side), show.legend = F)
    # colored labels for each side
    gg <- gg + geom_text(data = baselabels, aes(x = x, y = y, hjust = just, vjust = 0,
                                                label = text, color = side),
                         family = font, size = label_size, show.legend = F)
    # branch labels: GO terms
    gg <- gg + geom_text(data = treedata, aes(x = xst + text_x_offset,
                                              y = yst + (text_x_offset * tanpi(angle)) + text_y_offset,
                                              label = text,angle = textangle * 180,
                                              hjust = just), family = font, vjust = 0, size = label_size, show.legend = F)
    # p-value axis at the bottom
    gg <- gg + geom_segment(data = baseaxis, aes(x = xst, y = yst, xend = xst + len * cospi(angle),
                                                 yend = yst + len * sinpi(angle)),
                            size = 0.4, color = "black", show.legend = F)
    # labels for p-value axis
    if (plot_pvalue_labels) {
        gg <- gg + geom_text(data = baseaxis, aes(x = xst, y = yst - 1.2,
                                                              hjust = 1, vjust = 1, label = text),
                             family=font, size = label_size, angle = 45, show.legend = F, color = "black", parse = T)
    }
    # root box at the bottom
    gg <- gg + geom_label(data = base, aes(x = xst, y = yst, label = text),
                          size = label_size, nudge_y = 0.1, vjust = 1, parse = b_parse,color = root_color, family=font)
    # scale_y_continuous expands the plot in y direction; 10% at the bottom and 20% at the top
    gg <- gg + scale_size_identity() +
      scale_color_manual(values = colors, guide = F, drop = F) +
      scale_fill_manual(values = colors, guide = F, drop = F) +
      scale_y_continuous(expand = expand_scale(mult = c(0.2, 0.2))) +
      scale_x_continuous(expand = expand_scale(mult = expand_x)) +
      coord_fixed(ratio = 1) + theme_void()
    ggsave(filename = out_file, plot = gg, width = plot_size, height = plot_size)

}
