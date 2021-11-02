;;; fftw.lisp

;;; A simple FFTW wrapper, based on code from Patrick Klein.

(in-package :approxfun.coremath)

(eval-when
    (:compile-toplevel :load-toplevel :execute)

  (cffi:define-foreign-library libfftw3
    #+:x86_64 (:darwin "libfftw3_64.dylib")
    #-:x86_64 (:darwin "libfftw3.dylib")
    (:unix (:or "libfftw3.so.3" "libfftw3.so"))
    (t (:default "libfftw3")))

  (cffi:use-foreign-library libfftw3))

(cffi:defcenum fftw-direction
  (:forward -1)
  (:backward 1))

(cffi:defcenum fftw-planner-flags
  (:measure  #x000000)
  (:destroy-input #x000001)
  (:unaligned #x000002)
  (:converse-memory #x000004)
  (:exhaustive #x000008)
  (:preserve-input #x000010)
  (:patient  #x000020)
  (:estimate #x000040)
  (:estimate-patient #x000080)
  (:believe-pcost #x000100)
  (:no-dft-r2hc #x000200)
  (:no-nonthreads #x000400)
  (:no-buffering #x000800)
  (:no-indirect-op #x001000)
  (:allow-large-generic #x002000)
  (:no-rank-splits #x004000)
  (:no-vrank-splits #x008000)
  (:no-vrecurse #x010000)
  (:no-simd #x020000)
  (:no-slow #x040000)
  (:no-fixed-radix-large-n #x080000)
  (:allow-pruning #x100000)
  (:wisdom-only #x200000))


(cffi:defcfun ("fftw_plan_dft_1d" fftw-plan-dft-1d :library libfftw3) :pointer
  (n :int)
  (in (:pointer :double))
  (out (:pointer :double))
  (sign :int)
  (flags :unsigned-int))

(cffi:defcfun ("fftw_execute" fftw-execute) :void
  (plan :pointer))

(cffi:defcfun ("fftw_destroy_plan" fftw-destroy-plan) :void
  (plan :pointer))

(defstruct plan
  c-in
  c-out
  c-plan)

(defun create-plan (n sign flags)
  "Construct a FFTW execution plan for 1D DFT of length N."
  (let ((sign (cffi:foreign-enum-value 'fftw-direction sign))
        (flags (reduce #'logior flags
                       :initial-value 0
                       :key #'(lambda (ff)
                                (cffi:foreign-enum-value 'fftw-planner-flags
                                                         ff)))))
    (let ((plan (make-plan :c-in (cffi:foreign-alloc :double
                                                     :count (* 2 n))
                           :c-out (cffi:foreign-alloc :double
                                                      :count (* 2 n)))))
      (setf (plan-c-plan plan)
            (fftw-plan-dft-1d n (plan-c-in plan) (plan-c-out plan)
                              sign flags))
      plan)))

(defun execute (plan in out)
  "Execute a PLAN, using complex arrays IN and OUT for storage."
  (declare (type plan plan)
           (type (simple-array (complex double-float) (*)) in out)
           (optimize (speed 3)))
  (let ((c-in (plan-c-in plan)))
    (loop :for kk :of-type fixnum :from 0 :below (length in)
          :for rr :of-type fixnum :from 0 :by 2
          :for ii :of-type fixnum :from 1 :by 2
          :do (setf (cffi:mem-aref c-in :double rr) (realpart (aref in kk))
                    (cffi:mem-aref c-in :double ii) (imagpart (aref in kk)))))
  (fftw-execute (plan-c-plan plan))
  (let ((c-out (plan-c-out plan)))
    (loop :for kk :of-type fixnum :from 0 :below (length out)
          :for rr :of-type fixnum :from 0 :by 2
          :for ii :of-type fixnum :from 1 :by 2
          :do (let ((cc (complex (cffi:mem-aref c-out :double rr)
                                 (cffi:mem-aref c-out :double ii))))
                (setf (aref out kk) cc)))))

(defun destroy-plan (plan)
  "Release resources associated with PLAN."
  (fftw-destroy-plan (plan-c-plan plan))
  (cffi:foreign-free (plan-c-out plan))
  (cffi:foreign-free (plan-c-in plan)))


(defun fft (data &key (direction :forward))
  "Compute the Fast Fourier Transform of DATA.

Keyword Arguments:
  DIRECTION - either :FORWARD or :BACKWARDS, indicating the FFT or its inverse.
"
  (declare (type (simple-array (complex double-float) (*)) data))
  (let ((plan (create-plan (length data) direction '(:estimate)))
        (out (make-array (length data) :element-type '(complex double-float))))
    (unwind-protect 
         (execute plan data out)
      (destroy-plan plan))
    out))
