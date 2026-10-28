(require [amb3 [*]])
(import amb3 [amb3 [*]]
        [xarray :as xr]
        [matplotlib :as mpl]
        math pickle)

;;; netcdf containing horizontal grid data

(defn xr-open [fname]
  (xr.open-dataset fname :decode-cf False :decode-times False))

(defn str-ctypes [s] (ctypes.c-char-p (bytes s :encoding "ascii")))
(defn char-ctypes [s] (ctypes.c-char (bytes s :encoding "ascii")))

(defclass Lib []
  (defn --init-- [me]
    (sv me.lib (npy.ctypeslib.load-library "libmpas" ".")))
  (defn unit-test [me]
    (me.lib.unit-test))
  (defn set-mesh [me d]
    (sv c-int ctypes.c-int arr npy.array asc as-ctypes
        me.c-on-c (arr d.cellsOnEdge) me.e-on-c (arr d.edgesOnCell)
        me.v-on-c (arr d.verticesOnCell) me.c-on-e (arr d.cellsOnEdge)
        me.ne-on-c (arr d.nEdgesOnCell) me.v-on-e (arr d.verticesOnEdge)
        me.c-on-v (arr d.cellsOnVertex) me.e-on-v (arr d.edgesOnVertex)
        me.e-on-e (arr d.edgesOnEdge) me.ne-on-e (arr d.nEdgesOnEdge)
        me.w-on-e (arr d.weightsOnEdge)
        me.xedg (arr d.xEdge) me.yedg (arr d.yEdge) me.zedg (arr d.zEdge)
        me.xcel (arr d.xCell) me.ycel (arr d.yCell) me.zcel (arr d.zCell)
        me.xvtx (arr d.xVertex) me.yvtx (arr d.yVertex) me.zvtx (arr d.zVertex)
        me.nedge (first me.v-on-e.shape)
        (, me.ncell me.max-e-on-c) me.e-on-c.shape
        me.max-e-on-e (second me.w-on-e.shape)
        me.nvtx (len me.xvtx) me.nvdegree (second me.e-on-v.shape)
        me.dvEdge (arr d.dvEdge) me.dcEdge (arr d.dcEdge)
        me.angleEdge (arr d.angleEdge) me.areaCell (arr d.areaCell))
    (me.lib.set-mesh
     (c-int me.ncell) (c-int me.nedge) (c-int me.nvtx) (c-int me.max-e-on-c)
     (c-int me.nvdegree)
     (asc me.ne-on-c) (asc me.c-on-c) (asc me.e-on-c) (asc me.v-on-c)
     (asc me.c-on-e) (asc me.v-on-e)
     (asc me.c-on-v) (asc me.e-on-v)
     (asc me.xcel) (asc me.ycel) (asc me.zcel)
     (asc me.xedg) (asc me.yedg) (asc me.zedg)
     (asc me.xvtx) (asc me.yvtx) (asc me.zvtx)
     (asc me.dvEdge) (asc me.dcEdge)
     (asc me.angleEdge) (asc me.areaCell)
     ;; for edge velocity reconstruction
     (asc me.e-on-e) (asc me.ne-on-e) (asc me.w-on-e) (c-int me.max-e-on-e)))
  (defn check [me d &optional [time 0] [lvl 0] [out 42]]
    (me.set-mesh d)
    (sv c-int ctypes.c-int asc as-ctypes
        out-arr (npy.zeros me.nedge))
    (me.lib.check (c-int out) (asc out-arr))
    out-arr)
  (defn make-graph-file [me d outname]
    (me.set-mesh d)
    (me.lib.make-graph-file (str-ctypes outname)))
  (defn get-oriented-cells [me d]
    (me.set-mesh d)
    (sv oriented-elems (npy.zeros me.v-on-c.shape :dtype ctypes.c-int))
    (me.lib.get-oriented-cells (as-ctypes oriented-elems))
    oriented-elems)
  (defn get-oriented-cell-vertex-lonlat-list [me d &optional wrap]
    (svifn wrap True)
    (me.set-mesh d)
    (sv ll-list (npy.zeros (* 2 me.max-e-on-c me.ncell)))
    (me.lib.get-oriented-cell-vertex-lonlat-list (as-ctypes ll-list)
                                                 (ctypes.c-int wrap))
    ll-list)
  (defn get-oriented-cell-vertex-xy-list [me d]
    (me.set-mesh d)
    (sv ll-list (npy.zeros (* 2 me.max-e-on-c me.ncell))
        dbl ctypes.c-double)
    (me.lib.get-oriented-cell-vertex-xy-list (dbl d.x-period) (dbl d.y-period)
                                             (as-ctypes ll-list))
    ll-list)
  (defn make-yz-image [me d nz f3d]
    (me.set-mesh d)
    (sv nya (npy.zeros 1 :dtype npy.int32)
        f2d (npy.zeros f3d.size))
    (me.lib.make-yz-image (ctypes.c-int nz) (as-ctypes f3d)
                          (as-ctypes nya) (as-ctypes f2d))
    (sv ny (first nya))
    (.transpose (.reshape (get f2d (slice 0 (* ny nz))) (, ny nz))))
  (defn vert-slice-get-cols [me d cgids &optional axis]
    (assert (or (none? axis) (= (len axis) 3)))
    (me.set-mesh d)
    (sv coords-list (npy.zeros (* 2 (len cgids)))
        a (ctypes.c-int (len cgids))
        b (as-ctypes (npy.array cgids :dtype npy.int32))
        c (as-ctypes coords-list))
    (if (none? axis)
      (me.lib.vert-lat-slice-get-cols a b c)
      (me.lib.vert-axis-slice-get-cols a (as-ctypes axis) b c))
    (npy.reshape coords-list (, (len cgids) 2)))
  (defn test-siqp [me d]
    (me.set-mesh d)
    (me.lib.test-siqp))
  (defn test-ir [me d &optional [output False]]
    (me.set-mesh d)
    (if output
        (sv data (npy.zeros (* 4 me.ncell)))
        (sv data (npy.zeros 1)))
    (me.lib.unittest-ir (ctypes.c-int output) (as-ctypes data))
    data)
  (defn ir [me d &optional [output False] [optstr "ng"] [nstep 120]]
    (me.set-mesh d)
    (if output
      (sv ntracer (if (= (second optstr) "a") 4 1)
          n (* ntracer me.ncell)
          ini (npy.zeros n) mid (npy.zeros n) fin (npy.zeros n))
      (sv ini (npy.zeros 1) mid (npy.zeros 1) fin (npy.zeros 1)))
    (me.lib.advect-ir (ctypes.c-int output)
                      (str-ctypes optstr) (ctypes.c-int nstep)
                      (as-ctypes ini) (as-ctypes mid) (as-ctypes fin))
    (, ini mid fin))
  (defn ir-movie [me d &optional [optstr "ng"] [nstep 120]]
    (assert (!= (second optstr) "a"))
    (me.set-mesh d)
    (sv output 2 ini (npy.zeros (* me.ncell (inc nstep)))
        mid (npy.zeros 1) fin (npy.zeros 1))
    (me.lib.advect-ir (ctypes.c-int output)
                      (str-ctypes optstr) (ctypes.c-int nstep)
                      (as-ctypes ini) (as-ctypes mid) (as-ctypes fin))
    ini)
  (defn advect-mesh [me d &optional [optstr "na2200"] [nstep 8] [stepno 1]]
    (me.set-mesh d)
    (sv edge-ctrs (= (get optstr 2) "3")
        extra-factor (if edge-ctrs 2 1)
        n (* 2 extra-factor me.max-e-on-c me.ncell)
        ll-orig (npy.zeros n) ll-advected (npy.zeros n))
    (me.lib.advect-mesh (str-ctypes optstr) (ctypes.c-int nstep) (ctypes.c-int stepno)
                        (as-ctypes ll-orig) (as-ctypes ll-advected))
    (, ll-orig ll-advected (* 2 extra-factor me.max-e-on-c)))
  (defn integrate-gallery-velocity-field [me vf-short lons0 lats0 zetas0 ntime]
    (sv nspace (len lats0)
        ts (npy.zeros (* nspace ntime))
        m (fn [x]
            (vectorize (npy.dot (col-vec (npy.array x))
                                (row-vec (npy.ones ntime)))))
        lats (m lats0) lons (m lons0) zetas (m zetas0)
        r (fn [x] (.reshape x (, nspace ntime))))
    (me.lib.integrate-gallery-velocity-field
      (char-ctypes vf-short) (ctypes.c-int nspace) (ctypes.c-int ntime)
      (as-ctypes ts) (as-ctypes lons) (as-ctypes lats) (as-ctypes zetas))
    (, (r ts) (r lons) (r lats) (r zetas))))

(defn poly-prep-z [z &optional clim]
  (sv z (vectorize (npy.array z))
      znnan (get z (npy.logical-not (npy.isnan z)))
      zmin (min znnan) zmax (max znnan)
      zmin (if (npy.isnan zmin) 0 zmin)
      zmax (if (npy.isnan zmax) zmin zmax)
      zmax (if (= zmax zmin) (+ zmin 1e-3) zmax)
      zmin (if (none? clim) zmin (first clim))
      zmax (if (none? clim) zmax (last clim)))
  (, z zmin zmax))

(defn poly-prep-colors [z zmin zmax &optional cmap-obj]
  (svifn cmap-obj (pl.get-cmap "viridis"))
  (sv cmap (if (hasattr cmap-obj "colors")
             (npy.array cmap-obj.colors)
             (cmap-obj (range 0 cmap-obj.N)))
      nc (first cmap.shape)
      cmap-idxs (npy.array (npy.floor (* (/ (- z zmin)
                                            (- zmax zmin)) nc))
                           :dtype int)
      (get cmap-idxs (npy.isnan z)) 0
      (get cmap-idxs (<  cmap-idxs 0 )) 0
      (get cmap-idxs (>= cmap-idxs nc)) (dec nc)
      clrs (get cmap cmap-idxs))
  (, cmap-obj clrs))

(defn poly-draw [ax vtxs z &optional show-mesh colorbar clim cmap-obj]
  (svifn show-mesh False colorbar True)
  (pl.axes ax)
  (sv (, z zmin zmax) (poly-prep-z z clim)
      (, cmap-obj clrs) (poly-prep-colors z zmin zmax :cmap-obj cmap-obj)
      p (matplotlib.collections.PolyCollection vtxs))
  ;;(.set-array p (npy.array z))
  ;;(.set-facecolor p clrs)
  (.set-color p clrs)
  (when show-mesh (.set-edgecolor p "k"))
  (.add-collection ax p)
  (when colorbar
    (pl.colorbar (mpl.cm.ScalarMappable (mpl.colors.Normalize zmin zmax)
                                        cmap-obj)))
  (, zmin zmax))

(defclass PolyPlotter []
  (sv latlon 0 xy 1)

  (defn --init-- [me d &optional dims]
    (svifn dims me.latlon)
    (me.%init-c++-impl d dims))

  (defn %init-hy-impl [me d dims]
    (assert (= dims me.latlon))
    (sv u (Lib)
        elems (u.get-oriented-cells d)
        x d.lonVertex y d.latVertex
        nelem (first elems.shape)
        me.vtxs [])
    (for [ie (range nelem)]
      (sv ne (get u.ne-on-c ie)
          rows (slice 0 ne)
          elem (get (get elems ie) rows)
          xy (npy.zeros (, ne 2))
          lon (get x elem)
          f (* 0.25 math.pi)
          split (> (- (max lon) (min lon)) f))
      (when split
        (sv mask (< lon (- f math.pi)))
        (+= (get lon mask) (* 2 math.pi)))
      (sv (get xy (, rows 0)) lon
          (get xy (, rows 1)) (get y elem))
      (.append me.vtxs xy)))

  (defn %init-c++-impl [me d dims]
    (sv u (Lib)
        ll-list (case/eq dims
                         [me.latlon (u.get-oriented-cell-vertex-lonlat-list d)]
                         [me.xy (u.get-oriented-cell-vertex-xy-list d)]
                         [:else (assert False)])
        me.vtxs [])
    (for [ie (range u.ncell)]
      (sv ne (get u.ne-on-c ie)
          xy (npy.zeros (, ne 2))
          osx (* ie 2 u.max-e-on-c)
          osy (+ osx u.max-e-on-c)
          rows (slice 0 ne)
          (get xy (, rows 0)) (get ll-list (slice osx (+ osx ne)))
          (get xy (, rows 1)) (get ll-list (slice osy (+ osy ne))))
      (.append me.vtxs xy)))

  (defn draw [me ax z &optional show-mesh colorbar clim cmap-obj]
    (poly-draw ax me.vtxs z show-mesh colorbar clim cmap-obj)))

(defclass VertPlotter []
  (defn --init-- [me d cgids nvert &optional axis]
    (me.%init-c++-impl d cgids nvert axis))
  
  (defn %init-c++-impl [me d cgids nvert axis]
    (sv u (Lib)
        me.thetas-list (u.vert-slice-get-cols d cgids axis)
        ncell (len cgids)
        me.vtxs [])
    (for [ie (range ncell)]
      (sv (, hlo hhi) (get me.thetas-list ie))
      (for [k (range nvert)]
        (sv vlo (/ k nvert) vhi (/ (inc k) nvert))
        (.append me.vtxs (npy.array [[hlo vlo] [hhi vlo]
                                     [hhi vhi] [hlo vhi]])))))
  
  (defn draw [me ax z &optional show-mesh colorbar clim]
    (poly-draw ax me.vtxs z show-mesh colorbar clim)))

(defn draw-mesh [u extra-vtxs ll]
  (sv vtxs-fac (if extra-vtxs 2 1)
      n-per-h (* vtxs-fac u.max-e-on-c)
      n-per (* 2 n-per-h)
      npoly (// (len ll) n-per))
  (assert (zero? (% (len ll) n-per)))
  (for [i (range npoly)]
    (sv ne (* vtxs-fac (get u.ne-on-c i))
        os (* i n-per)
        xs (get ll (slice os (+ os ne)))
        os (+ os n-per-h)
        ys (get ll (slice os (+ os ne))))
    (pl.plot xs ys "k-")
    (pl.plot [(last xs) (first xs)] [(last ys) (first ys)] "k-")))

;;; general

(defn xyz->latlon [x y z]
  (sv f (fn [e] (** (first e) 2))
      r (math.sqrt (+ (f x) (f y) (f z)))
      q (/ z r)
      mask (> q 1) (get q mask) 1
      mask (< q -1) (get q mask) -1)
  (, (npy.arcsin q) (npy.arctan2 y x)))

(defn format-time [totsec]
  (assert (= (int totsec) totsec))
  (sv totsec (int totsec)
      hour (math.floor (/ totsec 3600))
      min (% (math.floor (/ totsec 60)) 60)
      sec (% totsec 60))
  (.format "{:02d}:{:02d}:{:02d}" hour min sec))

(expect (format-time 172233) "47:50:33")

;;; standalone mpas

(defn get-analytical [res-tag time-slice]
  (import analytical)
  (analytical.tracer "/compyfs/brad323/mpas/ocean/cosine_bell"
                     res-tag time-slice))

(defn analyze-parse [fname &optional [d None]]
  (svifn d {})
  (sv txt (.split (readall fname) "\n"))
  (for [ln txt]
    (when (< (len ln) 4) (break))
    (sv pre (cut ln 0 10)
        toks (.split ln))
    (cond [(in ">>" pre)
           (sv optstr (second toks))]
          [(in "re l2" pre)
           (sv idx (int (get toks 3))
               l2err (float (get toks 4))
               optstr (+ (first optstr) (get (, "t" "g" "c" "s") idx) (cut optstr 2)))]
          [(in "re mass" pre)
           (sv masserr (float (get toks 2)))]
          [(in "lim" pre)
           (sv toks toks
               limerr (max (float (get toks 1)) (float (get toks 2))))
           (assoc-nested-append d (, optstr 'l2err) l2err)
           (assoc-nested-append d (, optstr 'masserr) masserr)
           (assoc-nested-append d (, optstr 'limerr) limerr)]))
  d)

(defn analyze-replace [d drep]
  (for [k (.keys drep)]
    (assoc d k (get drep k)))
  d)

(defn analyze-plot [d outfn &optional [bare False]]
  (sv seqs (.keys d))
  (with [(pl-plot (, 9 6) outfn)]
    (for [sp (range 3)]
      (pl.subplot 1 3 (inc sp))
      (sv flow (get "rdn" sp))
      (for [s seqs]
        (if-not (= (first s) flow) (continue))
        (when (in (cut s 4) (, "l0" "0c" "lc"))
          (sv limerr (get d s 'limerr) masserr (get d s 'masserr))
          (when (> (max limerr) 1e-14) (prf " limerr {} {}" s limerr))
          (when (> (max masserr) 1e-12) (prf "masserr {} {}" s masserr)))
        (dont unless (in (cut s 4) (, "00" "0c" "lc")) (continue))
        (sv pat (+ (get {"t2" "k" "g2" "b" "c2" "g" "s2" "r"
                         "t3" "r" "g3" "c" "c3" "y" "s3" "k"}
                        (cut s 1 3))
                   (get {"1" "s" "2" "o" "t" "v"} (get s 3))
                   (get {"00" ":" "l0" "-." "0c" "--" "lc" "-"} (cut s 4)))
            y (get d s 'l2err))
        (pl.plot (npy.array (list (range 1 (inc (len y)))))
                 (npy.log2 (npy.array y))
                 pat :fillstyle "none"
                 :label (if (= (get s 1) "t") (cut s 2))))
      (pl.yticks (range -24 0))
      (axis-tight-pad)
      (pl.legend :loc "best")
      (my-grid))))

(defn study-plot [d outfn &optional [bare False]]
  (sv seqs (.keys d))
  (with [(pl-plot (, 9 6) outfn)]
    (for [sp (range 3)]
      (pl.subplot 1 3 (inc sp))
      (sv flow (get "rdn" sp))
      (for [s seqs]
        (if-not (= (first s) flow) (continue))
        (when (in (cut s 4) (, "l0" "0c" "lc"))
          (sv limerr (get d s 'limerr) masserr (get d s 'masserr))
          (when (> (max limerr) 1e-14) (prf " limerr {} {}" s limerr))
          (when (> (max masserr) 1e-12) (prf "masserr {} {}" s masserr)))
        (sv s1 s)
        (when (= (len s) 5) (sv s1 (+ s "0")))
        (sv pat (+ (get {"t" "k" "g" "r" "c" "b" "s" "g"} (get s 1))
                   ;;(get {1 "s" 2 "v" 3 "o"} (int (get s 2)))
                   (if (> (len s) 3)
                       (get {"1" "s" "2" "v" "t" "*" "b" "o"} (get s 3))
                       ".")
                   (if (> (len s) 3)
                       (get {"e00" "-" "200" "--" "e0c" "-." "20c" ":"
                             "10c" "--" "t0c" "-" "b0c" "-."
                             "100" "--" "t00" "-" "b00" "-."} (cut s1 3))
                       ""))
            y (get d s 'l2err))
        (pl.plot (npy.array (list (range 1 (inc (len y)))))
                 (npy.log2 (npy.array y))
                 pat :fillstyle "none"
                 :label (if (= (get s 1) "t") (cut s 2))))
      (pl.yticks (range -24 0))
      (axis-tight-pad)
      (pl.legend :loc "best")
      (my-grid))))

(defn nice-plot [d outfn &optional [bare False] [format None]]
  (svifn format "pdf")
  (sv seqs (.keys d))
  (with [(pl-plot (, 12 6) outfn :format format)]
    (for [sp (range 3)]
      (pl.subplot 1 3 (inc sp))
      (sv flow (get "rdn" sp))
      (for [s seqs]
        (if-not (= (first s) flow) (continue))
        (sv pat (+ (get {"t" "k" "g" "b" "c" "g" "s" "r"} (get s 1))
                   (get {"00" "s:" "l0" "x-." "0c" ".--" "lc" "o-"} (cut s 4)))
            y (get d s 'l2err))
        (pl.semilogy (npy.array (list (range 1 (inc (len y)))))
                     (npy.array y)
                     pat :fillstyle "none"
                     :label (cond [(and (zero? sp) (= (get s 1) "t"))
                                   (get {"00" "none" "l0" "limiter" "0c" "cedr" "lc" "both"}
                                        (cut s 4))]
                                  [(and (= sp 1) (= (cut s 4) "lc"))
                                   (get {"t" "XYZTrig" "g" "GH" "c" "CB" "s" "SC"}
                                        (get s 1))])))
      (pl.xticks [1 2 3 4 5] ["480km\n48 steps" "240km\n96 steps" "120km\n192 steps"
                              "60km\n384 steps" "30km\n768 steps"])
      ;;(axis-tight-pad)
      (my-grid)
      (pl.title (get (, "Rotation" "Divergent flow" "Nondivergent flow") sp))
      (when (zero? sp)
        (pl.ylabel "$l_2$ error"))
      (when (< sp 2) (pl.legend :loc "best")))))

;; graph file (mpas .graph.info) for gpmetis
(defn make-graph-file [ncname outname]
  (with [d (xr-open ncname)]
    (sv u (Lib))
    (u.make-graph-file d outname)))

(if-main

(sv amb3.*when-inp-verbosity* 2)

(sv *outdir* "/qfs/people/brad323/v2/www/random/mpas/")

(when-inp ["unittest"]
  (sv u (Lib))
  (u.unit-test))

(when-inp ["tracer1" {:res str :time int}]
  (with [d (xr.open-dataset
            (.format (+ "/compyfs/brad323/mpas/ocean/cosine_bell/"
                        "{}/default/forward/output/output.0001-01-01_00.00.00.nc")
                     res))]
    (sv lat d.latCell lon d.lonCell
        t1 (get d.tracer1 (, time (s-all) 0)))
    (print (len lat) (len lon) (len t1))
    (sv analytical (get-analytical res time))
    (with [(pl-plot (, 6 9)
                    (+ *outdir* (.format "t1-{}-{}" res time))
                    :format "png")]
      (for [i (range 3)]
        (pl.subplot 3 1 (inc i))
        (pl.tripcolor lon lat
                      (case/eq i
                               [0 t1]
                               [1 analytical]
                               [2 (- t1 analytical)])
                      :shading "gouraud")
        (pl.axis "equal")
        (case/eq 1
                 [0 (pl.xlim (, 2.5 4)) (pl.ylim (, -.4 .4))]
                 [1 (pl.xlim (, 2.4 3.9)) (pl.ylim (, -.4 .4))])
        (pl.colorbar)))))

(when-inp ["check-ll" {:fname str}]
  (sv arr npy.array)
  (defn prerr [tlat tlon lat lon]
    (sv tlat (arr tlat) tlon (arr tlon)
        laterr (npy.max (npy.abs (- lat tlat)))
        lonerr 0)
    (for [i (range (len lon))]
      (sv a (get lon i) b (- (get tlon i)) twopi (* 2 math.pi)
          lonerr (max lonerr
                      (min (abs (+ a b))
                           (abs (+ a twopi b))
                           (abs (+ a (- twopi) b))))))
    (prf "laterr {:1.3e} lonerr {:1.3e}" laterr lonerr))
  (with [d (xr.open-dataset fname)]
    (sv ncell (first d.edgesOnCell.shape)
        nedge (first d.verticesOnEdge.shape)
        (, clat clon) (xyz->latlon (arr d.xCell) (arr d.yCell) (arr d.zCell))
        (, elat elon) (xyz->latlon (arr d.xEdge) (arr d.yEdge) (arr d.zEdge)))
    (prerr d.latCell d.lonCell clat clon)
    (prerr d.latEdge d.lonEdge elat elon)))

(when-inp ["plot-areaCell" {:fname str}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        lat d.latCell lon d.lonCell)
    (with [(pl-plot (, 18 9) "fig/areaCell" "png")]
      (pl.tripcolor lon lat (/ d.areaCell (npy.max d.areaCell)) :shading "flat")
      (pl.axis "equal")
      (pl.colorbar))))

(when-inp ["plot-cc" {:fname str}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        clat d.latCell clon d.lonCell)
    (with [(pl-plot (, 18 9) "fig/cc" "png")]
      (pl.plot clon clat "k.")
      (pl.axis "equal"))))

(when-inp ["plot-nedge" {:fname str}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        clat d.latCell clon d.lonCell nedge d.nEdgesOnCell)
    (with [(pl-plot (, 18 9) "fig/nedge" "png")]
      (pl.tripcolor clon clat nedge :shading "flat")
      (pl.axis "equal"))))

(when-inp ["plot-edge-angle" {:fname str}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        err (u.check d :out 0)
        lat d.latEdge lon d.lonEdge)
    (with [(pl-plot (, 18 9) "fig/angleEdge-err" "png")]
      (pl.tripcolor lon lat (npy.log10 err) :shading "gouraud")
      (pl.axis "equal")
      (pl.colorbar))))

(when-inp ["check" {:fname str}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib))
    (u.check d)))

(when-inp ["check-ad-hoc" {:fname str}]
  ;; want to run some q&d checks on mpaso.WC14to60E2r3.rstFromG-anvil.201002.nc
  (with [d (xr.open-dataset fname)]
    (sv npa npy.array
        ncell (npy.size d.latCell))
    (prf "ncell {}" ncell)
    (prf "min cellsOnEdge {}" (npy.min (npa d.cellsOnEdge)))
    (prf "max cellsOnEdge {}" (npy.max (npa d.cellsOnEdge)))
    (prf "number of 0 entries in cellsOnEdge {}"
         (npy.sum (npa (zero? (npa d.cellsOnEdge)) :dtype int)))
    (prf "min areaCell {}" (npy.min (npa d.areaCell)))
    (prf "max areaCell {}" (npy.max (npa d.areaCell)))))

(sv sset (fn [a b c] {:separate-rho a :supercycle b :nstepfac c})
    *sc-settings* (, (sset 0 0 1)   ; baseline
                     (sset 0 0 3)   ; cmp with subsequent
                     (sset 1 1 3)   ; rho but no supcyc
                     (sset 1 3 1))  ; supcyc
    *mesh-dir* "/compyfs/brad323/mpas/meshes/"
    *mesh-name-template* "mesh.QU.{:d}km.151026.nc"
    *km-seq* (, 480 240 120 60)
    *nstep-480* {"n" 150 "d" 90 "r" 75})

(defn make-settings-label [s]
  (.format "rho {} supcyc {} nstepfac {}"
           (:separate-rho s) (:supercycle s) (:nstepfac s)))

(defn save-lims [] (, (pl.xlim) (pl.ylim)))
(defn restore-lims [lims] (pl.xlim (first lims)) (pl.ylim (second lims)))

(when-inp ["ir" {:fname str :optstr str :nstep int :output int}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        patches True
        (, ic mid fin) (u.ir d :output output :optstr optstr :nstep nstep))
    (when output
      (if patches
          (sv pp (PolyPlotter d))
          (sv lon d.lonCell lat d.latCell))
      (with [(pl-plot (, 18 9) "fig/ir-test" "png")]
        (for [i (range 4)]
          (sv ax (pl.subplot 2 2 (inc i))
              data (case/eq i [0 ic] [1 mid] [2 fin] [3 (- fin ic)]
                            [-3 (npy.log10 (npy.abs (- fin ic)))]))
          (if patches
              (pp.draw ax data)
              (do
                (pl.tripcolor lon lat data :shading "gouraud")
                (pl.colorbar)))
          (dont when (< i 3)
            (case/eq (second optstr) ["t" (pl.clim 0 1)] ["g" (pl.clim 0 0.95)]
                     ["c" (pl.clim 0.05 1)] ["s" (pl.clim 0.1 1)]))
          (pl.axis "tight"))))))

(when-inp ["advect-mesh" {:fname str :optstr str :nstep int :stepno int}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib) pi math.pi
        extra-vtxs (= (get optstr 2) "3")
        (, llo lla n-per) (u.advect-mesh d :optstr optstr :nstep nstep :stepno stepno))
    (with [(pl-plot (, 9 9) "fig/advect-mesh" "png")]
      (for [i (range 2)]
        (sv ax (pl.subplot 2 1 (inc i)))
        (draw-mesh u extra-vtxs (case/eq i [0 llo] [1 lla]))
        (pl.xlim (- pi) pi)
        (pl.ylim (* -0.5 pi) (* 0.5 pi))))))

(when-inp ["ir-test" {:fname str :output int}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        data (u.test-ir d :output output)))
  (when output
    (sv lon d.lonCell lat d.latCell)
      (with [(pl-plot (, 18 9) "fig/ir-test" "png")]
        (for [i (range 4)]
          (pl.subplot 2 2 (inc i))
          (pl.tripcolor lon lat
                        (cut data (* u.ncell i) (* u.ncell (inc i)))
                        :shading "gouraud")
          (pl.axis "tight")
          (pl.colorbar)))))

(when-inp ["ir-movie" {:fname str :optstr str :nstep int}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        patches True
        frames (u.ir-movie d :optstr optstr :nstep nstep))
    (if patches
        (sv pp (PolyPlotter d))
        (sv lon d.lonCell lat d.latCell))
    (for [step (range (inc nstep))]
      (with [(pl-plot (, 9 4.5) (.format "fig/ir-movie-{:03d}" step) "png")]
        (sv ax (pl.subplot 1 1 1)
            data (cut frames (* step u.ncell) (* (inc step) u.ncell)))
        (if patches
            (pp.draw ax data)
            (do
              (pl.tripcolor lon lat data :shading "gouraud")
              (pl.colorbar)))
        (pl.text -3 1.3
                 (get {"r" "Solid-body rotation"
                       "n" "Nondivergent flow"
                       "d" "Divergent flow"}
                      (first optstr))
                 :color "w" :fontsize 15)
        (pl.axis "equal")
        (pl.axis "off")))))

(when-inp ["siqp-test" {:fname str}]
  (with [d (xr.open-dataset fname)]
    (sv u (Lib))
    (u.test-siqp d)))

(when-inp ["cmp-trajectories" {:fname str :optstr str :nstep int :output int}]
  (defn reldif [w a b]
    (math.sqrt (npy.sum (* w
                           (/ (npy.sum (** (- (npy.array b) (npy.array a)) 2))
                              (npy.sum (** (npy.array a) 2)))))))
  (with [d (xr.open-dataset fname)]
    (sv u (Lib)
        (, eic emid efin) (u.ir d :output output :optstr (+ optstr "e") :nstep nstep)
        (, tic tmid tfin) (u.ir d :output output :optstr (+ optstr "2") :nstep nstep))
    (when output
      (prf "re l2 mid {:1.3e} end {:1.3e}"
           (reldif d.areaCell tmid emid) (reldif d.areaCell tfin efin))
      (sv lon d.lonCell lat d.latCell)
      (with [(pl-plot (, 18 9) "fig/ir-test" "png")]
        (for [i (range 4)]
          (pl.subplot 2 2 (inc i))
          (pl.tripcolor lon lat
                        (case/eq i [0 (- tic eic)] [1 (- tmid emid)] [2 (- tfin efin)]
                                 [3 (- (- tfin tic) (- efin eic))])
                        :shading "gouraud")
          (pl.axis "tight")
          (pl.colorbar))))))

(when-inp ["dev-ice"]
  (sv m (get-python-mod "vs.py" "tmp")
      vs (npy.transpose (npy.array m.vs)))
  (with [(pl-plot (, 6 6) "fig/vs")]
    (cond [False
           (sv ax ((. (pl.figure) gca) :projection "3d"))
           (ax.view-init 90 270)
           (ax.plot (get vs 0) (get vs 1) (get vs 2) "k-")
           (pl.xlabel "x") (pl.ylabel "y")
           (dont pl.axis "off")]
          [True
           (pl.plot (get vs 0) (get vs 1) "k-")
           (pl.axis "equal")])))

(when-inp ["study-o3"]
  (sv d (analyze-parse "data/study-o3-7.txt" :d None))
  (do
    (sv keys (list (.keys d)))
    (dont for [k keys]
      (unless (in (cut k 2) (, "220c" "3e0c" "3e00"))
        (del (get d k)))))
  (study-plot d "fig/study-o3" :bare False))

(when-inp ["study-time"]
  (sv d (analyze-parse "data/study-time-8.txt" :d None)
      d1 (analyze-parse "data/study-time-9.txt" :d None)
      d (analyze-replace d d1))
  (sv d (analyze-parse "data/study-time-9c.txt" :d None))
  (study-plot d "fig/study-time-9c" :bare False))

(when-inp ["analyze"]
  (sv sfx "-9"
      d (analyze-parse (.format "data/runall{}.txt" sfx))
      )
  (dont
    (sv keys (list (.keys d)))
    (for [k keys]
      (unless (in (cut k 2) (, "220c" "320c"))
        (del (get d k)))))
  (analyze-plot d (.format "fig/runall{}" sfx) :bare False)
  ;;(return)
  (nice-plot d (.format "fig/nice{}" sfx) :bare False :format "pdf"))

(defn ir-analyze-parse [fname]
  (sv txt (.split (readall fname) "\n")
      d {})
  (for [ln txt]
    (sv toks (.split ln " "))
    (cond [(in ">>" ln)
           (sv mesh (last toks)
               mesh-type (if (in ".QU." mesh) 'qu 'jigsaw)
               km -1)
           (for [i (, 240 120 60)]
             (when (in (.format "{:d}" i) mesh) (sv km i) (break)))]
          [(in "optstr" ln)
           (sv optstr (get toks 1))]
          [(in "re l2" ln)
           (sv l2-err (float (last toks)))
           (assoc-nested d (, optstr mesh-type km) l2-err)]))
  d)

(defn ir-analyze-plot [d-ir fname]
  (sv odes "rn"
      reses (, 240 120 60)
      x [1 2 3])
  (with [(pl-plot (, 9 12) fname :format "png")]
    (for [(, di d) (enumerate (, d-ir))]
      (sv sp-os (if (zero? di) 0 2))
      (for [(, oi ode) (enumerate odes)]
        (pl.subplot 2 2 (+ sp-os (inc oi)))
        (for [(, optstr d1) (.items d)]
          (unless (= (first optstr) ode) (continue))
          (for [mesh-type (, 'qu 'jigsaw)]
            (sv y [])
            (for [res reses]
              (.append y (get d1 mesh-type res)))
            (pl.plot x (npy.log2 (npy.array y))
                     (+ (get {"t" "k" "g" "r" "c" "b"} (get optstr 1))
                        (get {"2" "v" "3" "o" "4" "s"} (get optstr 2))
                        (if (= mesh-type 'qu) "-" "--"))
                     :label (+ optstr "-" (get {'qu "qu" 'jigsaw "un"} mesh-type))
                     :fillstyle "none")))
        (pl.xticks x reses)
        (pl.yticks (npy.array (list (range -17 0))))
        (my-grid)
        (axis-tight-pad)
        (pl.legend :loc "best" :ncol 2)
        (when (zero? oi) (pl.ylabel "$\log_2$ $l_2$ error"))
        (pl.title (get {"r" "Rotation" "n" "Nondivergent flow"} ode))
        (when (zero? oi)
          (sv yl (pl.ylim))
          (pl.text 0.5 (+ (second yl) .7)
                   (if (zero? di)
                       "Eulerian flux-form methods"
                       "RF-CISL method")
                   :fontsize 14
                   :horizontalalignment "left"))))))

(when-inp ["ir-analyze"]
  (sv d-ir (ir-analyze-parse "data/run-ir-ir-equiv-0.txt"))
  (ir-analyze-plot d-ir "fig/ir-run-0"))

;; hy mpas.hy make-graph-file $mesh.nc $mesh.graph.info
;; gpmetis $mesh.graph.info $nparts:int
(when-inp ["make-graph-file" {:ncname str :outname str}]
  (make-graph-file ncname outname))

) ;; if-main
