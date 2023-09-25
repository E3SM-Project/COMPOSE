(require [amb3 [*]])
(import amb3 [amb3 [*]]
        poly
        [scipy.linalg :as linalg])

(defn make-gll-mesh [xb np &optional [end False]]
  (sv ne (dec (len xb))
      d (dec np)
      n (* ne d)
      xgll (npy.zeros (if end (inc n) n))
      ref (* 0.5 (+ 1 (npy.array (cut (poly.get-gll-x np) 0 -1)))))
  (for [ie (range ne)]
    (sv dx (- (get xb (inc ie)) (get xb ie))
        (get xgll (slice (* ie d) (* (inc ie) d))) (+ (get xb ie) (* dx ref))))
  (when end (sv (get xgll n) (+ (- (last xb) (first xb)) (first xgll))))
  xgll)

(defn cond-2norm [A]
  (sv s (linalg.svd A :compute-uv False))
  (/ (max s) (min s)))

(defn np-pad [a arr b]
  (npy.append a (npy.append arr b)))

(defn make-cells [ncell &optional [pattern 'uniform] [perturb-factor 0.5]]
  (defn dx->xb [dx]
    (sv xb (npy.cumsum dx)
        xb (npy.append 0 xb)
        xb (/ xb (get xb -1)))
    xb)
  (defn make-sin [perturb-factor]
    (sv dx (npy.sin (* 1 math.pi (npy.linspace 0 1 ncell)))
        dx (+ 2 dx)
        dx (* (+ 1 (* perturb-factor (npy.random.rand ncell))) dx))
    (dx->xb dx))
  (assert (> ncell 1))
  (case/eq pattern
           ['perturb
            (sv r (npy.random.rand (dec ncell))
                xb (+ (npy.linspace 0 1 (inc ncell))
                      (np-pad 0 (* perturb-factor (- r 0.5) (/ ncell)) 0)))]
           ['sin (sv xb (make-sin perturb-factor))]
           ['random
            (sv xb (np-pad 0 (npy.random.rand (dec ncell)) 1))
            (.sort xb)]
           ['random-refine
            (sv xb (make-cells (// ncell 2) 'random)
                xb (npy.concatenate (, xb (cb->cc xb))))
            (.sort xb)]
           ['uniform
            (sv xb (npy.linspace 0 1 (inc ncell)))]
           ['uniformbut1
            (sv xb (npy.linspace 0 1 (inc ncell))
                (get xb 1) (+ (get xb 1) (/ 0.5 ncell)))]
           ['alternating
            (sv nper 3
                dx (npy.zeros (, nper (// (+ ncell nper -1) nper)))
                (get dx (, 0 (s-all))) 0.75
                (get dx (, 1 (s-all))) 0.65
                (get dx (, 2 (s-all))) 0.6
                dx (cut (.reshape (.transpose dx) (npy.prod dx.shape)) 0 ncell)
                xb (dx->xb dx))]
           ['alternating-perturb
            (sv nper 3
                dx (npy.zeros (, nper (// (+ ncell nper -1) nper))))
            (for [i (range nper)]
              (sv (get dx (, i (s-all)))
                  (+ 1 (* 2 (npy.random.rand) perturb-factor))))
            (sv dx (cut (.reshape (.transpose dx) (npy.prod dx.shape)) 0 ncell)
                xb (dx->xb dx))]
           [:else (raisefmt "Not a pattern: {:s}" (name pattern))])
  (when False
    (sv dx-max (npy.max (npy.diff xb))
        dx-min (npy.min (npy.diff xb)))
    (print dx-max dx-min (/ dx-max dx-min)))
  (assert (npy.all (> (cut xb 1) (cut xb 0 -1))))
  xb)

(defn make-ic [x shape-list &optional [width 0.2]]
  (in-require (np-array? x))

  (when (in 'locality shape-list)
    (when (> (len shape-list) 1)
      (raisefmt "make-ic: locality can be used only on its own."))
    (sv width 0.2
        hw (/ width 2)))

  (sv bot 0.1
      top 1.0
      y (* bot (npy.ones (len x)))
      nshapes (len shape-list)
      hs (/ (* 2 nshapes))
      hw (/ width 2)
      x0 (* 1 hs))
  (defn draw-trig [ctr]
    (sv (get y (s-all)) (+ 1 (npy.sin (* 2 math.pi x)))))
  (defn make-mask [ctr]
    (* (>= x (- ctr hw)) (<= x (+ ctr hw))))
  (defn draw-delta [ctr]
    (sv mask (>= x ctr)
        i (first (get (npy.array (range (len x))) mask))
        (get y i) top))
  (defn draw-rect [ctr]
    (sv (get y (make-mask ctr)) (* 0.8 top)))
  (defn draw-tri [ctr]
    (sv left  (* (>= x (- ctr hw)) (<= x ctr))
        right (* (<= x (+ ctr hw)) (>= x ctr))
        lx (get x left) rx (get x right)
        l0 (- ctr hw) r0 ctr
        slope (/ (- top bot) hw)
        (get y left)  (+ bot (* slope (- lx l0)))
        (get y right) (+ top (* slope (- r0 rx)))))
  (defn draw-sin [ctr]
    (sv mask (make-mask ctr)
        xm (get x mask)
        cosxm (npy.cos (* math.pi (/ (- xm ctr) hw)))
        (get y mask) (+ bot
                        (* 0.5 (- (* 0.95 top) bot)
                           (+ 1 cosxm)))))
  (defn draw-gauss [ctr]
    (sv (cut y 0)
        (+ y
           (* (- top bot)
              (npy.exp (- (/ (** (- x ctr) 2)
                            (** (* 0.27 width) 2))))))))
  (defn draw-locality [ctr]
    (sv mask (make-mask ctr)
        not-mask (= mask False)
        (get y mask) (* 1 top)
        (get y not-mask) (* 0.25
                            (+ 1 (npy.sin
                                  (* 2 math.pi (get x not-mask)))))))
  ;; treat the line as a circle embedded in a plane from which a 2D gaussian is
  ;; sampled. this gives a truly C^inf periodic function with substantial
  ;; bandwidth (the second unlike trig).
  (defn draw-gauss2d [ctr]
    (sv theta (* 2 math.pi x)
        theta-ctr (* 2 math.pi ctr)
        u (npy.cos theta) v (npy.sin theta)
        u0 (npy.cos theta-ctr) v0 (npy.sin theta-ctr))
    (sv (cut y 0)
        (+ y
           (* (- top bot)
              (npy.exp (- (/ (+ (** (- u u0) 2) (** (- v v0) 2))
                             (** (* 2 width) 2))))))))

  (for [(, i shape) (enumerate shape-list)]
    ((case/eq shape
              ['delta draw-delta]
              ['rect draw-rect]
              ['tri draw-tri]
              ['sin draw-sin]
              ['gauss draw-gauss]
              ['gauss2d draw-gauss2d]
              ['trig draw-trig]
              ['locality
               (sv width 0.2 hw (* 0.5 width))
               draw-locality]
              [:else (raise (Exception "Not a shape."))])
     (+ x0 (* i 2 hs))))
  y)

(defn cc [x]
  (* 0.5 (+ (cut x 0 -1) (cut x 1))))

(defn duplicate-dg-data-as-cg [data ne &optional [end False]]
  (sv np (len data)
      n (* (dec np) ne)
      a (npy.zeros (if end (inc n) n)))
  (for [ie (range ne)]
    (sv (get a (slice (* ie (dec np)) (* (inc ie) (dec np)))) (cut data 0 -1)))
  (when end
    (sv (get a -1) (last data)))
  a)

(defn get-domain-of-dependence [A]
  (sv n (get A.shape 1)
      dod [])
  (for [row A]
    (sv nz [])
    (for [j (range n)]
      (if (!= (get row j) 0)
        (.append nz j)))
    (.sort nz)
    (.append dod nz))
  dod)

(defn frame-init []
  (pl.show :block False))
(defn frame-flip [U &optional [clim None]]
  (pl.clf)
  (pl.imshow U :interpolation "none")
  (unless (none? clim)
    (pl.clim clim))
  (pl.pause 0.05))

(defn caas [dod w yp y]
  (sv n (len yp)
      l (npy.zeros n)
      u (npy.zeros n))
  (for [i (range n)]
    (sv yp-i (get yp (get dod i))
        lo (min yp-i)
        (get l i) lo
        hi (max yp-i)
        (get u i) hi
        (get y i) (min hi (max lo (get y i)))))
  (sv mass-prev (sum (* w yp))
      mass-curr (sum (* w y))
      dm (- mass-prev mass-curr))
  (if (>= dm 0)
    (sv capacity (- u y)
        alpha (/ dm (sum (* w capacity)))
        y (+ y (* alpha capacity)))
    (sv capacity (- l y)
        alpha (/ dm (sum (* w capacity)))
        y (+ y (* alpha capacity))))
  y)
