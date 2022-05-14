(in-package :approxfun)

(defmethod interactive-gnuplot::translate-to-fragment ((obj interval))
  (fragment "[~F:~F]" (interval-lower obj) (interval-upper obj)))

(defun plot (apfun &key
                     (num-samples 500)
                     (points nil)
                     (xrange nil)
                     (yrange nil))
  "Plot APFUN.

Keyword Arguments:
- NUM-SAMPLES is the number of function evaluations to use
- POINTS is a list of additional points (cons x y) to plot
- XRANGE is an interval indicating the range for the x axis
- YRANGE is an interval indicating the range for the y axis"
  (check-type apfun chebyshev-approximant)
  (let* ((fn (chebyshev-approximant-interp-fn apfun))
         (d (chebyshev-approximant-interval apfun))
         (l (interval-lower d))
         (r (interval-upper d))
	 (y-min most-positive-double-float)
	 (y-max most-negative-double-float))
    (uiop:with-temporary-file (:stream data :pathname tmp)
      (loop :for x :from l :to r :by (cl:/ (- r l) num-samples)
            :for y := (funcall fn x)
            :do (format data "~&~F ~F" x y)
		(setf y-min (min y y-min)
		      y-max (max y y-max)))
      (finish-output data)
      (uiop:with-temporary-file (:stream data-points :pathname tmp-pts)
        (loop :for (x . y) :in (coerce points 'list)
              :do (format data-points "~&~F ~F" x y)
		  (setf y-min (min y y-min)
			y-max (max y y-max)))
        (finish-output data-points)
	(interactive-gnuplot:gnuplot
	    (:set :background :rgb "white")
	    (:set :key :off)
	    (:set :xrange (or xrange d))
	    (:set :yrange (or yrange
			      (let ((fudge-factor (* 0.1 (- y-max y-min))))
				(interval (- y-min fudge-factor) (+ y-max fudge-factor)))))
	    (:plot (fragment "'~A' w lines~@[, '~A' w points~]"
			     tmp (and points tmp-pts))))))))

(defun plot-coefficients (apfun)
  "Plots the absolute value of the Chebyshev coefficients of APFUN to FILE.

If DISPLAY is T, then this will open the resulting plot."
  (check-type apfun chebyshev-approximant)
  (let ((coeffs (chebyshev-approximant-coeffs apfun))
	(y-max 0d0))
    (uiop:with-temporary-file (:stream data :pathname tmp)
      (loop :for i :from 0 :below (length coeffs)
            :for c := (abs (aref coeffs i))
            :do (format data "~&~D ~F" i c)
		(setf y-max (max y-max c)))
      (finish-output data)
      (gnuplot
	(:set :background :rgb "white")
	(:set :key :off)
	(:set :yrange (interval 1d-16 y-max))
	(:set :xrange (interval 0 (length coeffs)))
	(:set :logscale :y)
	(:plot (fragment "'~A'" tmp))
	(:unset :logscale :y)))))
