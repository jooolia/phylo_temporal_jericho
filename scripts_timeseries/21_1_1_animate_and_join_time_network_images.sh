### make the plots over time into little animationss

convert -adjoin `ls -1v ../figures/overall_no_delay_t_*.png` ../figures/overall_no_delay_time_plot.pdf

## so figure does not keep getting split up every time I run the make file. 
touch ../figures/overall_no_delay_time_plot.png

montage -mode concatenate -density 300 `ls -1v ../figures/overall_no_delay_t_*.png` ../figures/overall_no_delay_time_plot_montage.png

## the -lv puts them in the correct over
convert -delay 70 `ls -1v ../figures/overall_no_delay_t_*.png` ../figures/overall_no_delay_time.gif

convert ../figures/overall_no_delay_time.gif -delay 70 -layers TrimBounds -coalesce -gravity south \
          -pointsize 24 -fill '#FFF8' -stroke '#0008' \
          \( -clone 0 -annotate 0 "2010-06-22" \) -swap 0 +delete \
          \( -clone 1 -annotate 0 "2010-07-06" \) -swap 1 +delete \
          \( -clone 2 -annotate 0 "2010-07-20" \) -swap 2 +delete \
          \( -clone 3 -annotate 0 "2010-08-05" \) -swap 3 +delete \
          \( -clone 4 -annotate 0 "2010-08-17" \) -swap 4 +delete \
          \( -clone 5 -annotate 0 "2010-08-31" \) -swap 5 +delete \
          \( -clone 6 -annotate 0 "2010-09-14"  \) -swap 6 +delete \
          \( -clone 7 -annotate 0 "2010-09-28" \) -swap 7 +delete \
          \( -clone 8 -annotate 0 "2010-10-12" \) -swap 8 +delete \
          \( -clone 9 -annotate 0 "2010-10-26"  \) -swap 9 +delete \
          \( -clone 10 -annotate 0 "2010-11-09" \) -swap 10 +delete \
          \( -clone 11 -annotate 0 "2010-11-23" \) -swap 11 +delete \
          \( -clone 12 -annotate 0 "2010-12-06"  \) -swap 12 +delete \
          \( -clone 13 -annotate 0 "2010-12-20" \) -swap 13 +delete \
          \( -clone 14 -annotate 0 "2011-01-10" \) -swap 14 +delete \
          \( -clone 15 -annotate 0 "2011-01-18" \) -swap 15 +delete \
          \( -clone 16 -annotate 0 "2011-01-29"  \) -swap 16 +delete \
          \( -clone 17 -annotate 0 "2011-01-31" \) -swap 17 +delete \
          \( -clone 18 -annotate 0 "2011-02-02" \) -swap 18 +delete \
          \( -clone 19 -annotate 0 "2011-02-04" \) -swap 19 +delete \
          \( -clone 20 -annotate 0 "2011-02-06" \) -swap 20 +delete \
          \( -clone 21 -annotate 0 "2011-02-08" \) -swap 21 +delete \
          \( -clone 22 -annotate 0 "2011-02-10"  \) -swap 22 +delete \
          \( -clone 23 -annotate 0 "2011-02-24" \) -swap 23 +delete \
          \( -clone 24 -annotate 0 "2011-03-10" \) -swap 24 +delete \
          \( -clone 25 -annotate 0 "2011-03-24" \) -swap 25 +delete \
          \( -clone 26 -annotate 0 "2011-04-07" \) -swap 26 +delete \
          \( -clone 27 -annotate 0 "2011-04-21" \) -swap 27 +delete \
          \( -clone 28 -annotate 0 "2011-05-05"  \) -swap 28 +delete \
          \( -clone 29 -annotate 0 "2011-05-18" \) -swap 29 +delete \
          \( -clone 30 -annotate 0 "2011-06-07"  \) -swap 30 +delete \
          \( -clone 31 -annotate 0 "2011-06-21" \) -swap 31 +delete \
          \( -clone 32 -annotate 0 "2011-07-05" \) -swap 32 +delete \
          \( -clone 33 -annotate 0 "2011-07-19" \) -swap 33 +delete \
          -layers OptimizeFrame   ../figures/no_delay_animated_with_date_.gif


#-fill '#FFF8' removed this to make it easier to see. 
convert -density 300 `ls -1v ../figures/overall_no_delay_t_*.png` -coalesce -gravity south \
          -pointsize 24  -stroke '#0008' \
          \( -clone 0 -annotate 0 "2010-06-22" \) -swap 0 +delete \
          \( -clone 1 -annotate 0 "2010-07-06" \) -swap 1 +delete \
          \( -clone 2 -annotate 0 "2010-07-20" \) -swap 2 +delete \
          \( -clone 3 -annotate 0 "2010-08-05" \) -swap 3 +delete \
          \( -clone 4 -annotate 0 "2010-08-17" \) -swap 4 +delete \
          \( -clone 5 -annotate 0 "2010-08-31" \) -swap 5 +delete \
          \( -clone 6 -annotate 0 "2010-09-14"  \) -swap 6 +delete \
          \( -clone 7 -annotate 0 "2010-09-28" \) -swap 7 +delete \
          \( -clone 8 -annotate 0 "2010-10-12" \) -swap 8 +delete \
          \( -clone 9 -annotate 0 "2010-10-26"  \) -swap 9 +delete \
          \( -clone 10 -annotate 0 "2010-11-09" \) -swap 10 +delete \
          \( -clone 11 -annotate 0 "2010-11-23" \) -swap 11 +delete \
          \( -clone 12 -annotate 0 "2010-12-06"  \) -swap 12 +delete \
          \( -clone 13 -annotate 0 "2010-12-20" \) -swap 13 +delete \
          \( -clone 14 -annotate 0 "2011-01-10" \) -swap 14 +delete \
          \( -clone 15 -annotate 0 "2011-01-18" \) -swap 15 +delete \
          \( -clone 16 -annotate 0 "2011-01-29"  \) -swap 16 +delete \
          \( -clone 17 -annotate 0 "2011-01-31" \) -swap 17 +delete \
          \( -clone 18 -annotate 0 "2011-02-02" \) -swap 18 +delete \
          \( -clone 19 -annotate 0 "2011-02-04" \) -swap 19 +delete \
          \( -clone 20 -annotate 0 "2011-02-06" \) -swap 20 +delete \
          \( -clone 21 -annotate 0 "2011-02-08" \) -swap 21 +delete \
          \( -clone 22 -annotate 0 "2011-02-10"  \) -swap 22 +delete \
          \( -clone 23 -annotate 0 "2011-02-24" \) -swap 23 +delete \
          \( -clone 24 -annotate 0 "2011-03-10" \) -swap 24 +delete \
          \( -clone 25 -annotate 0 "2011-03-24" \) -swap 25 +delete \
          \( -clone 26 -annotate 0 "2011-04-07" \) -swap 26 +delete \
          \( -clone 27 -annotate 0 "2011-04-21" \) -swap 27 +delete \
          \( -clone 28 -annotate 0 "2011-05-05"  \) -swap 28 +delete \
          \( -clone 29 -annotate 0 "2011-05-18" \) -swap 29 +delete \
          \( -clone 30 -annotate 0 "2011-06-07"  \) -swap 30 +delete \
          \( -clone 31 -annotate 0 "2011-06-21" \) -swap 31 +delete \
          \( -clone 32 -annotate 0 "2011-07-05" \) -swap 32 +delete \
          \( -clone 33 -annotate 0 "2011-07-19" \) -swap 33 +delete \
          -layers OptimizeFrame   ../figures/no_delay_montage_with_date.png

montage -mode concatenate -density 300 `ls -1v ../figures/no_delay_montage_with_date-*.png` ../figures/no_delay_time_plot_montage_labelled.png