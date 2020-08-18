(in-package :approxfun)

(defparameter *display-command-name* "open"
  "The command to use when displaying a plot.")

(defun plot (apfun file &key
                      (display nil)
                      (num-samples 500))
  "Plot APFUN to FILE.

Uses NUM-SAMPLES function evaluations. If DISPLAY is T, then this will open the
resulting plot.
"
  (check-type apfun chebyshev-approximant)
  (let ((fn (chebyshev-approximant-interp-fn apfun))
	(d (chebyshev-approximant-domain apfun)))
    (uiop:with-temporary-file (:stream data :pathname tmp)
      (loop :for x :from (domain-lower d) :to (domain-upper d) :by (/ 2 num-samples)
            :do (format data "~&~F ~F" x (funcall fn x)))
      (finish-output data)
      (uiop:run-program
       (list "gnuplot" "-e"
             (format nil "set terminal svg; set key off; set xrange [~,3F:~,3F]; plot '~A' w lines"
		     (domain-lower d) (domain-upper d)
                     tmp))
       :output file)))
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
             (format nil "set terminal svg; set key off; set logscale y; plot '~A"
                     tmp))
       :output file)))
  (when display
    (uiop:run-program (list *display-command-name* (uiop:unix-namestring (truename file)))))
  file)
