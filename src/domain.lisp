(in-package #:approxfun)

(define-condition domain-error (simple-error)
  ()
  (:documentation "An error associated with an invalid domain."))

(defun domain-error (format-control &rest format-args)
  "Signal a DOMAIN-ERROR, with the associated message."
  (error 'domain-error
         :format-control format-control
         :format-arguments format-args))

(defun make-domain (lower upper)
  (let ((interval (interval lower upper)))
    (values interval
            (approxfun #'identity :interval interval))))

(defgeneric domain (obj)
  (:documentation "Get the domain associated to OBJ."))

(defun domain= (da db)
  (interval= da db))
