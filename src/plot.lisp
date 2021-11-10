(in-package :approxfun)

(defparameter *display-command-name* "open"
  "The command to use when displaying a plot.")

(defparameter *plot-debug-stream* nil
  "An optional stream to dump plotter commands to.")

(defun plot (apfun file &key
                          (display nil)
                          (num-samples 500)
                          (points nil)
                          (xrange nil)
                          (yrange nil))
  "Plot APFUN to FILE.

Keyword Arguments:
- If DISPLAY is T, then the resulting plot will be opened.
- NUM-SAMPLES is the number of function evaluations to use
- POINTS is a list of additional points (cons x y) to plot
- XRANGE is an interval indicating the range for the x axis
- YRANGE is an interval indicating the range for the y axis
"
  (declare (ignore yrange))
  (check-type apfun chebyshev-approximant)
  (let* ((fn (chebyshev-approximant-interp-fn apfun))
         (d (chebyshev-approximant-interval apfun))
         (l (interval-lower d))
         (r (interval-upper d)))
    (uiop:with-temporary-file (:stream data :pathname tmp)
      (loop :for x :from l :to r :by (cl:/ (- r l) num-samples)
            :for y := (funcall fn x)
            :do (format data "~&~F ~F" x y)
                (when *plot-debug-stream*
                  (format *plot-debug-stream* "~&~F ~F" x y)))
      (finish-output data)
      (uiop:with-temporary-file (:stream data-points :pathname tmp-pts)
        (loop :for (x . y) :in (coerce points 'list)
              :do (format data-points "~&~F ~F" x y))
        (finish-output data-points)
        (let ((gnuplot-command
                (with-output-to-string (gp)
                  (format gp "set terminal svg enhanced background rgb 'white'; ")
                  (format gp "set key off; ")
                  (format gp "set xrange [~,3F:~,3F]; "
                          (if xrange (interval-lower xrange) (interval-lower d))
                          (if xrange (interval-upper xrange) (interval-upper d)))
                  (format gp "plot '~A' w lines~@[, '~A' w points~]"
                          tmp (and points tmp-pts)))))
          (when *plot-debug-stream*
            (format *plot-debug-stream* "~%~A~%" gnuplot-command))
          (uiop:run-program
           (list "gnuplot" "-e" gnuplot-command)
           :output file)))))
  (when display
    (uiop:run-program (list *display-command-name* (uiop:unix-namestring (truename file)))))
  file)

(defun plot-coefficients (apfun file &key (display nil))
  "Plots the absolute value of the Chebyshev coefficients of APFUN to FILE.

If DISPLAY is T, then this will open the resulting plot."
  (check-type apfun chebyshev-approximant)
  (let ((coeffs (chebyshev-approximant-coeffs apfun)))
    (uiop:with-temporary-file (:stream data :pathname tmp)
      (loop :for i :from 0 :below (length coeffs)
            :for c := (aref coeffs i)
            :do (format data "~&~D ~F" i (abs c)))
      (finish-output data)
      (uiop:run-program
       (list "gnuplot" "-e"
             (format nil "set terminal svg enhanced background rgb 'white'; set key off; set logscale y; plot '~A"
                     tmp))
       :output file)))
  (when display
    (uiop:run-program (list *display-command-name* (uiop:unix-namestring (truename file)))))
  file)
