utMA <- function(values,times,tau=24)
{
# utMA = unequal time-step (simple) moving average
# code modified from: http://www.eckner.com/papers/ts_alg.pdf
# update (2017-08-28) http://www.eckner.com/papers/Algorithms%20for%20Unevenly%20Spaced%20Time%20Series.pdf
# default tau assumes times are in hours and observations occur approximately 3/day

    # if class(times) == 'difftime' convert to numeric
    if(class(times) == 'difftime')    times <- as.numeric(times)

    out <- rep(0,length(values))
    left <- 1; 
    roll_area <- left_area <- values[1] * tau
    out[1] <- values[1]; 
    for (right in 2:length(values))
    {
        # Expand interval on right end
        roll_area <- roll_area + values[right-1] * (times[right] - times[right-1]);

        # Remove truncated area on left end 
        roll_area = roll_area - left_area;

        # Shrink interval on left end 
        t_left_new = times[right] - tau; 
        while (times[left] <= t_left_new) 
        {
            roll_area = roll_area - values[left] * (times[left+1] - times[left]);
            left = left + 1;
        }

        # Add truncated area on left end
        left_area <- values[max(1, left-1)] * (times[left] - t_left_new)
        
        roll_area <- roll_area + left_area;
            
        # Save SMA value for current time window
        out[right] <- roll_area / tau; 
    }
    
    # should figure out how to replace values where no averaging occurs with original value
    out
}


trapezoid <- function(x1, x2, x3, y1, y3) 
{
    if (x2 == x3 | x2 < x1)
    {
        return((x3 - x2) * y1)
    } else {
        weight <- (x3 - x2) / (x3 - x1);
        y2 <- y1 * weight + y3 * (1 - weight); 
        return((x3 - x2) * (y2 + y3) / 2);
    }
}

utMAlin <- function(values,times,tau=24)
{
# utMA = unequal time-step (simple) moving average
# code modified from: http://www.eckner.com/papers/ts_alg.pdf
# default tau assumes times are in hours and observations occur approximately 3/day

    # if class(times) == 'difftime' convert to numeric
    if(class(times) == 'difftime')    times <- as.numeric(times)

    out <- rep(0,length(values))
    out[1] <- values[1]

    left <- 1; 
    roll_area <- left_area <- values[1] * tau

    for (right in 2:length(values))
    {
        # Expand interval on right end
        roll_area <- roll_area + (values[right-1] + values[right])/2 * (times[right] - times[right-1]);

        # Remove truncated area on left end 
        roll_area = roll_area - left_area;

        # Shrink interval on left end 
        t_left_new = times[right] - tau; 
        while (times[left] <= t_left_new) 
        {
            roll_area = roll_area - (values[left] + values[left+1]) / 2 * (times[left+1] - times[left]);
            left = left + 1;
        }

        # Add truncated area on left end
        left_area <- trapezoid(times[max(1, left-1)], t_left_new, times[left],
            values[max(1, left-1)], values[left]); 
        
        roll_area <- roll_area + left_area;
            
        # Save SMA value for current time window
        out[right] <- roll_area / tau; 
    }
    
    # should figure out how to replace values where no averaging occurs with original value
    out
}

countToMA <- function(dat,col2avg='cell.count',tau=12,dig=3)
{
    # tau is time window used to cacluate average
    out <- dat
    out[,col2avg] <- round(utMAlin(dat[,col2avg],dat$time,tau=tau),dig)
    out
}

