(require [amb3 [*]])
(import amb3 [amb3 [*]]
        [figsutils [*]]
        math glob struct)

(assoc matplotlib.rcParams "savefig.dpi" 300)
(do (pl-require-type1-fonts))

;;; parse cmd, L, M, C slmmir output

(defn acc-parse [fname &optional map-nstepfac]
  (sv txt (.split (readall fname) "\n")
      bo None d {})
  (for [ln txt]
    (sv ln2 (cut ln 0 2))
    (cond [(and (= "cm" ln2) (in "cmd>" ln))
           (sv cmd ln c (parse-cmd ln :map-nstepfac map-nstepfac))]
          [(= ln2 "M ")
           (sv mp (parse-midpoint-check ln))
           (unless (none? (geton d #*(+ (cmd->key-base c)
                                        (, (:ic mp) cyc (:ne c) :M))))
             ;; handle repeated ic used for src-term OOA testing.
             (assoc mp :ic (+ (:ic mp) "-src")))
           (assoc-nested d (+ (cmd->key-base c) (, (:ic mp) cyc (:ne c) :M))
                         {:l1 (:l1 mp) :l2 (:l2 mp)})]
          [(= ln2 "C ")
           (cond [(in "cycle" ln) (sv (, - - cyc) (sscanf ln "s,s,i"))]
                 [(or (in "PASS" ln) (in "FAIL" ln))
                  (sv (, - pf) (sscanf ln "s,s"))
                  (when (and (= pf "FAIL") (!= (:mono c) "none") (<= cyc 10))
                    (prf "FAIL {}" cmd))]
                 [:else
                  (sv cl (parse-C ln))
                  (unless (none? (geton d #*(+ (cmd->key-base c)
                                               (, (:ic cl) cyc (:ne c) :C))))
                    ;; handle repeated ic used for src-term OOA testing.
                    (assoc cl :ic (+ (:ic cl) "-src")))
                  (assoc-nested d (+ (cmd->key-base c) (, (:ic cl) cyc (:ne c) :C))
                                cl)])]
          [(= ln2 "L ")
           (cond [(in "L file" ln)
                  (assoc-nested d (+ (cmd->key-base c)
                                     (, "cos" cyc (:ne c) :L :mixing-file))
                                (last (.split ln)))]
                 [(in "phimin" ln)
                  (sv (, - ic - l1 - l2 - li - phimin - phimax)
                      (sscanf ln "s,s,s,f,s,f,s,f,s,f,s,f"))
                  (assoc-nested d (+ (cmd->key-base c) (, ic cyc (:ne c) :Lerr))
                                {:l1 l1 :l2 l2 :li li :phimin phimin :phimax phimax})]
                 [:else
                  (sv bo (parse-bakeoff-diag bo ln (:timeint c)))
                  (when (:done bo)
                    (assoc-nested d (+ (cmd->key-base c) (, "cos" cyc (:ne c) :L)) bo)
                    (sv bo None))])]))
  d)

;;; print a long table of all accuracy results

(defn acc-print-txt-table [c d &optional]
  (sv dent (Box) dent.indentation 0 dent.delta 1)
  (defn indent [] (+= dent.indentation dent.delta))
  (defn dedent [] (-= dent.indentation dent.delta))
  (sv buf [])
  (defn msg [s] (.append buf (+ (* " " dent.indentation) s)))
  (defn msg-pop [] (unless (empty? buf) (.pop buf)))
  (defn msg-clear [] (.clear buf))
  (defn msg-dump []
    (for [e buf] (prf e))
    (msg-clear))
  (defmacro dent-fn [title &rest body]
    `(do (msg ~title)
         (indent)
         ~@body
         (msg-pop)
         (dedent)))

  (sv unstab-thr 0.9 unstab {})
  (defn study-unstab [e-all cdrglb cdrlcl np ic nstepfac ode timeint]
    (sv e None)
    (for [ne c.nes]
      (sv pe e e (geton e-all ne))
      (when (none? e) (continue))
      (unless (or (> (get e :C :l2) 0.9)
                  (and (not (none? pe)) (> (get e :C :l2) (get pe :C :l2))
                       (> (get e :C :l2) 1e-12)))
        (continue))
      (dont prf "u {:2d} {} {:2d} {:2d} {} {} {:7.1e}"
           np ic ne nstepfac cdrglb cdrlcl (get e :C :l2))
      (sv key (, np cdrglb cdrlcl)
          use True)
      (when (in key unstab)
        (sv pv (get unstab key))
        (unless (or (and (= (:ic pv) "slo") (!= ic "slo"))
                    (and (< ne (:ne pv))
                         (not (and (!= (:ic pv) "slo") (= ic "slo"))))
                    (and (not (= ic "slo"))
                         (> (get e :C :l2) (get pv :l2))
                         (< (get pv :l2) unstab-thr)))
          (sv use False)))
      (when (and (!= np 4)
                 (or (!= timeint "interp") (= cdrglb "none")
                     (!= ode "divergent") (!= ic "gau") (!= nstepfac 1)))
        (continue))
      (when use
        (assoc unstab key {:ic ic :ne ne :l2 (get e :C :l2)
                           :nstepfac nstepfac :ode ode}))))
  (defn print-unstab []
    (sv keys (list (.keys unstab)))
    (.sort keys :key first)
    (for [k keys]
      (sv e (get unstab k))
      (prf "u {:2d} {} {:<12s} {:2d} {:2d} {:<9s} {} {:7.1e}"
           (first k) (:ic e) (:ode e) (:ne e) (:nstepfac e) (second k)
           (last k) (:l2 e))))

  (defn pr [&rest args]
    (prfno (* " " dent.indentation))
    (print #*args))
  (defn praccv [pne ne vp v]
    (if (or (none? vp) (zero? vp) (zero? v))
        (prfno " {:7.1e} ( ----)" v)
        (prfno " {:7.1e} ({:5.2f})" v (first (ooa [vp v] :x [ne pne])))))
  (defn pracc [diagnostic norm pne pe ne e]
    (sv v (geton e diagnostic norm))
    (when (none? v) (return))
    (sv vp (unless (none? pe) (get pe diagnostic norm)))
    (praccv pne ne vp v))
  (defn print-table [e-all &optional [Ldiags False]]
    (sv pe None pne None)
    (for [ne c.nes]
      (sv e (geton e-all ne))
      (when (none? e) (continue))
      (when (< (get e :C :l2) 1e-13) (continue))
      (prfno "{}{:3d}" (* " " dent.indentation) ne)
      (pracc ':C ':masscons pne pe ne e)
      (for [diag (, :C :M)]
        (sv first True)
        (for [norm (, :l1 :l2 :li)]
          (when (and (= diag :M) (!= norm :l2)) (continue))
          (when (none? (geton e diag norm)) (continue))
          (when first
            (prfno " |")
            (sv first False))
          (pracc diag norm pne pe ne e)))
      (when Ldiags
        (sv bo (geton e :L))
        (unless (none? bo)
          (for [r (, 0 1)]
            (unless (or (zero? r) (in :me-mixing bo)) (continue))
            (prfno " |")
            (sv s (if (zero? r) :mixing :me-mixing)
                pmixing (unless (none? pe) (get pe :L s))
                mixing (get bo s))
            (for [k (, :lr :lu :lo)]
              (praccv pne ne (unless (none? pmixing) (get pmixing k))
                      (get mixing k))))))
      (when (> (get e :C :l2) 0.9) (prfno " UNSTABLE"))
      (print)
      (sv pe e pne ne)))

  (pr (+ "         ne    mass                l1              l2            linf"
         "            mid l2"))
  (for [method c.methods]
    (dent-fn method
      (for [timeint c.timeints]
        (dent-fn (.format "timeint {}" timeint)
          (for [ode c.odes]
            (dent-fn ode
              (for [cdrgl c.cdrs]
                (sv (, cdrglb cdrlcl) cdrgl)
                (dent-fn (.format "{} {}" cdrglb cdrlcl)
                  (for [nstepfac c.nstepfacs]
                    (dent-fn (.format "nstepfac {}" nstepfac)
                      (for [ic c.ics]
                        (dent-fn ic
                          (for [cyc c.cycs]
                            (dent-fn (.format "cycle {}" cyc)
                              (for [np c.nps]
                                (sv prefine (if (or (= np 4) (= timeint "exact"))
                                                0 5)
                                    e (geton d timeint ode nstepfac method
                                             cdrglb cdrlcl prefine np ic cyc))
                                (when (none? e) (continue))
                                (dent-fn (.format "np {}" np)
                                  (msg-dump)
                                  (pr method timeint ode cdrglb cdrlcl nstepfac
                                      ic cyc np)
                                  (print-table e :Ldiags (= ic "cos"))
                                  (study-unstab e cdrglb cdrlcl np ic nstepfac
                                                ode timeint)))))))))))))))))
  (print-unstab))

;;; accuracy figs

(defclass AccFig []
  (defn --init-- [me])

  (defn get-defaults [me &optional context]
    (svifn context (get-context))
    (sv c context)
    {:method "pcsl" :nstepfac 1 :timeint "interp" :ode "nondivergent" :ic "gau"
     :cdrglb "caas-node" :cdrlcl "caas" :prefine 5 :nonuni 0 :cyc 1 :measure :l2
     :nps c.nps :nes c.nes :pat-line "-" :pat-clr (.copy c.npclrs) :C-line :C
     :pat-mark (.copy c.npmarks) :fs 11 :lw 1 :markersize 4 :yform :log2
     :xticks :degrees :ooa-text False :filter-floor None :figsize (, 4 4)
     :ref-ooa-2 False :ref-cos-033 True :pg None})

  (defn plot [me ax d-all &optional [o None]]
    (svifn o (me.get-defaults))
    (sv npa npy.array
        d1 (geton d-all (:timeint o) (:ode o) (:nstepfac o) (:method o)
                  (:cdrglb o) (:cdrlcl o) (:prefine o))
        gray (* 0.5 (npy.ones 3))
        y-ext [1e10 -1e10])
    (when (none? d1) (return))
    (for [np (:nps o)]
      (sv keys (if (none? (:pg o))
                   (, np (:ic o) (:cyc o))
                   (if (fn? (:pg o))
                       (, ((:pg o) np) np (:ic o) (:cyc o))
                       (, (:pg o) np (:ic o) (:cyc o))))
          d2 (geton d1 #*keys))
      (when (none? d2) (continue))
      (sv x [] y [])
      (for [ne (:nes o)]
        (sv val (geton d2 ne (:C-line o) (:measure o)))
        (when (none? val) (continue))
        (unless (none? (:filter-floor o))
          (when (< val (:filter-floor o)) (continue)))
        (.append x ne)
        (.append y val)
        (sv (get y-ext 0) (min (get y-ext 0) val)
            (get y-ext 1) (max (get y-ext 1) val)))
      (svb (pat (+ (get (:pat-clr o) np) (get (:pat-mark o) np) (:pat-line o)))
           (xtform (npy.log2 (npa x)))
           ((, xticks xtick-strs)
            (case/eq (:xticks o)
                     [:degrees
                      (sv d (nes->degstrs x))
                      (, (npy.log2 (npa (:ne d))) (:deg d))]
                     [:else (, xtform x)]))
           ((, ytform y-lbl y-tfn)
            (case/eq (:yform o)
                     [:log2 (, (npy.log2 (npa y)) "$\log_2$" npy.log2)]
                     [:log10 (, (npy.log2 (npa y)) "$\log_{10}$" npy.log2)]
                     [:semilogy (, (npa y) "$\log_{10}$" identity)]))
           (pl-plot (if (= (:yform o) :semilogy) ax.semilogy ax.plot)))
      (pl-plot xtform ytform pat
               :lw (:lw o) :markersize (:markersize o) :fillstyle "none")
      (when (= np (first (:nps o)))
        (pl.xticks xticks xtick-strs :fontsize (:fs o))
        (pl.ylabel (+ y-lbl " " (get {:l1 "$l_1$" :l2 "$l_2$" :li "$l_{\infty}$"
                                      :masscons "Mass conservation"}
                                     (:measure o))
                      " relative error")
                   :fontsize (:fs o))
        (cond [(= (:yform o) :log2)
               (pl.yticks (list (range -40 10)) :fontsize (:fs o))])
        (when (and (:ref-cos-033 o) (= (:ic o) "cos") (= (:ode o) "nondivergent"))
          (pl-plot xtform (* 0.033 (npy.ones (len xtform))) "k-."
                   :zorder -10 :lw (:lw o) :color gray)))
      (when (:ooa-text o)
        (sv i (- (len x) 2))
        (pl.text (- (get xtform (inc i)) 0.1) (* (get ytform (inc i)) 2)
                 (.format "{:1.1f}"
                          (- (first (ooa (cut y i (+ i 2))
                                         :x (cut x i (+ i 2))))))))
      (when (= np (last (:nps o)))
        (when (= (:yform o) :semilogy)
          (sv (, y ys) (pl.yticks))
          (when (> (/ (last y) (first y)) 1e3)
            (for [i (range (len ys))]
              (sv (get ys i) (.format "{}" (int (math.log10 (get y i))))))
            (pl.yticks y ys)))
        (when (:ref-ooa-2 o)
          (sv ytref (y-tfn (* 0.7 (first y-ext) (** (/ (last x) (npa x)) 2))))
          (pl-plot xtform ytref "k:" :color gray))))
    y-ext)

  (defn legend [me ax entries &optional [o None] [nps-legend True] bbox]
    (svifn o (me.get-defaults))
    (sv xl (pl.xlim) yl (pl.ylim) hs [] delta 2)
    (unless (empty? entries)
      (for [e entries]
        (sv h (ax.plot (- (first xl) delta) (- (first yl) delta)
                       (first e) :label (second e) :fillstyle "none"
                       :lw (:lw o) :markersize (:markersize o)))
        (.extend hs h))
      (sv l1 (pl.legend :handles hs :fontsize (- (:fs o) 1)
                        :bbox-to-anchor (if (none? bbox) (, 0 0.08) bbox)
                        :loc "lower left" :framealpha 1))
      (ax.add-artist l1))
    (when nps-legend
      (sv hs [])
      (for [np (:nps o)]
        (sv h (ax.plot (- (first xl) delta) (- (first yl) delta)
                       (+ (get (:pat-clr o) np) (get (:pat-mark o) np))
                       :lw (:lw o) :markersize (:markersize o)
                       :label (.format "{}" np) :fillstyle "none"))
        (.extend hs h))
      (sv l2 (pl.legend :handles hs :fontsize (- (:fs o) 1)
                        :ncol (len (:nps o)) :bbox-to-anchor (, 0 -0.01 1 0)
                        :loc "lower left" :mode "expand" :framealpha 1))
      (ax.add-artist l2))
    (pl.xlim xl) (pl.ylim yl))

  (defn title [me s &optional [o None]]
    (svifn o (me.get-defaults))
    (pl.title s :fontsize (:fs o))))

(defn make-nps-string [nps]
  (sv s (.format "$n_p$ {}" (first nps)))
  (for [ne (cut nps 1)] (sv s (+ s (.format ", {}" ne))))
  s)

(defn nstepfac->word [nstepfac]
  (get {1 "long" 3 "medium" 5 "short"} nstepfac))

(defn make-title [main o &optional extra]
  (svifn extra "")
  (+ main "\n"
     (flow-short2long (:ode o)) ", "
     (ic-short2long (cut (:ic o) 0 3)) ", "
     (nstepfac->word (:nstepfac o)) " steps,\n"
     (if (= (:prefine o) 5) "$p$-refinement, " "")
     (+ (if (= (:cdrglb o) "none") "no " "") "property preservation")
     extra))

(defn fig-stab-cmp [c d]
  (sv p (AccFig)
      o (p.get-defaults c))
  (assoc o :ic "gau" :ode "divergent" :cdrglb "caas-node" :cdrlcl "caas"
         :measure :l2 :timeint "interp" :prefine 5 :yform :semilogy)
  (sv nps-str (make-nps-string (:nps o)))
  (with [(pl-plot (:figsize o) (+ c.fig-dir "stab-cmp-" (name (:measure o))))]
    (sv ax (pl.subplot 1 1 1))
    (assoc o :method "pcslu" :pat-line "-." :cyc 1) (p.plot ax d o)
    (assoc o :method "pcsl" :pat-line "--" :cyc 100) (p.plot ax d o)
    (assoc o :method "pcsl" :pat-line "-" :cyc 1 :ref-ooa-2 True) (p.plot ax d o)
    (if (= (:yform o) :semilogy)
        (pl.ylim (, 9e-7 1))
        (pl.ylim (, -20 0)))
    (my-grid)
    (p.legend ax (, (, "k-" "Islet 1 cycle") (, "k--" "Islet 100 cycles")
                    (, "k-." "Natural 1 cycle") (, "k:" "OOA 2")) :o o)
    (p.title (make-title "Islet stability:" o) o)))

(defn nextpow10 [f] (** 10.0 (math.ceil  (math.log10 f))))
(defn prevpow10 [f] (** 10.0 (math.floor (math.log10 f))))

(defn figs-acc [c d &optional prefix ref-ooa-2 legend
                general-timeint general-prefine show-linf pp]
  (svifn prefix "" ref-ooa-2 True legend True general-timeint "interp"
    general-prefine 5 show-linf True pp True)
  (sv p (AccFig)
      o (p.get-defaults c))
  (defn plot [o title plot-fn &optional [ref-ooa-2 False]]
    (sv fname (+ prefix "acc-" (:ode o) "-" (:ic o)
                 "-" (:timeint o) "-" (if (= (:cdrglb o) "none")
                                          "nopp" "pp")
                 "-fac" (str (:nstepfac o))))
    (print fname)
    (with [(pl-plot (:figsize o) (+ c.fig-dir fname))]
      (sv ax (pl.subplot 1 1 1))
      (plot-fn ax d o)
      (my-grid)
      (sv legs [(, "k-" "$l_2$") (, "k--" "$l_{\infty}$")])
      (when ref-ooa-2 (.append legs (, "k:" "OOA 2")))
      (when legend (p.legend ax legs :o o))
      (pl.ylabel "$\log_{10}$ relative error")
      (p.title (make-title title o) o)))
  (assoc o :ode "nondivergent" :ic "gau" :nstepfac 1
         :yform :semilogy :cdrglb "none" :cdrlcl "none"
         :timeint "exact" :prefine 0 :filter-floor 1e-11)
  (plot o "Islet empirical order of accuracy:"
        (fn [ax d o]
          (for [norm (, :l2 :li)]
            (assoc o :measure norm :pat-line (if (= norm :l2) "-" "--")
                   :ooa-text (= norm :li))
            (p.plot ax d o)
            (pl.ylim (, 1e-11 1))
            (sv e (npy.array [0 -2 -4 -6 -8 -10]))
            (pl.yticks (** 10.0 e) e))))
  (for [nstepfac (, 1 5)
        ic (, "gau" "cos" "slo")
        ode (, "nondivergent" "divergent" "rotate")]
    (assoc o :nstepfac nstepfac :ooa-text False :ode ode
           :cdrglb (if pp "caas-node" "none") :cdrlcl (if pp "caas" "none")
           :ic ic :filter-floor None :timeint general-timeint)
    (dont when (and (= ode "divergent") (= ic "gau")) (continue))
    (plot o "Islet accuracy:"
          (fn [ax d o]
            (sv ye [1e10 -1e10])
            (defn update-ye [ye1]
              (when (none? ye1) (return))
              (sv (get ye 0) (min (first ye) (first ye1))
                  (get ye 1) (max (last ye) (last ye1))))
            (for [norm (, :l2 :li)]
              (when (and (= norm :li) (not show-linf)) (continue))
              (assoc o :cdrglb (if pp "caas" "none") :timeint "exact" :prefine 0
                     :measure norm :pat-line (if (= norm :l2) "-" "--")
                     :ref-ooa-2 False)
              (sv ye1 (p.plot ax d o))
              (update-ye ye1)
              (assoc o :cdrglb (if pp "caas-node" "none")
                     :timeint general-timeint :prefine general-prefine
                     :ref-ooa-2 (and ref-ooa-2 (= norm :l2) (!= ic "slo")))
              (sv ye1 (p.plot ax d o))
              (update-ye ye1))
            (if (= ic "slo")
                (pl.ylim (, (* 0.7 (first ye)) (nextpow10 (second ye))))
                (pl.ylim (, (/ (first ye) 4) (nextpow10 (second ye))))))
          :ref-ooa-2 (and ref-ooa-2 (!= ic "slo")))))

(defn fig-acc-midpoint-check [c d]
  (sv p (AccFig)
      o (p.get-defaults c))
  (assoc o :ic "gau" :ode "nondivergent" :cdrglb "none" :cdrlcl "none"
         :measure :l2 :yform :semilogy :nstepfac 1)
  (sv nps-str (make-nps-string (:nps o)))
  (with [(pl-plot (:figsize o) (+ c.fig-dir "midpoint-check"))]
    (sv ax (pl.subplot 1 1 1))
    (assoc o :pat-line "-" :C-line :C :ooa-text True :prefine 0 :timeint "exact" :nps [4])
    (p.plot ax d o)
    (assoc o :prefine 5 :timeint "interp" :nps (cut c.nps 0 -1))
    (p.plot ax d o)
    (assoc o :prefine 5 :timeint "interp" :nps (cut c.nps -1) :ooa-text False)
    (p.plot ax d o)
    (assoc o :pat-line "--" :C-line :M :prefine 0 :timeint "exact" :nps [4])
    (p.plot ax d o)
    (assoc o :prefine 5 :timeint "interp" :nps c.nps)
    (p.plot ax d o)
    (pl.ylim (, 5e-10 1))
    (my-grid)
    (p.legend ax (, (, "k-" "1 cycle") (, "k--" "1/2 cycle")) :o o)
    (p.title (make-title "Trajectory interpolation:" o) o)))

(defn fig-acc-mimic-src-term-midpoint [c d &optional [np-minus-2 False]]
  (sv p (AccFig)
      o (p.get-defaults c))
  (defn plot [o title plot-fn &optional [ref-ooa-2 False]]
    (sv fname (+ "acc-pg-mimic-src-term-midpoint-" (:ode o) "-" (:ic o)
                 "-" (:timeint o) "-" (if (= (:cdrglb o) "none")
                                          "nopp" "pp")
                 "-fac" (str (:nstepfac o)) "-" (name (:measure o)))
        title (make-title title o :extra ", midpoint"))
    (print fname)
    (with [(pl-plot (:figsize o) (+ c.fig-dir fname))]
      (sv ax (pl.subplot 1 1 1))
      (plot-fn ax d o)
      (my-grid)
      (sv legs [(, "k-" "Reference") (, "k-." "$n_f=n_p$") (, "k--" "$n_f=2$")])
      (when np-minus-2 (.append legs (, "k:" "$n_f=n_p-2$")))
      (when ref-ooa-2 (.append legs (, "k:" "OOA 2")))
      (p.legend ax legs :o o)
      (pl.ylabel "$\log_{10}$ $l_2$ relative error")
      (p.title title o)))
  (for [nstepfac (, 5)
        ic (, "gau")
        ode (, "nondivergent")]
    (assoc o :nstepfac nstepfac :ooa-text False :ode ode :cdrglb "caas-node"
           :cdrlcl "caas" :ic ic :filter-floor None :norm :l2 :timeint "interp"
           :ref-cos-033 False :yform :semilogy :filter-floor 1e-11 :C-line :M)
    (plot o "Physics grid source term accuracy:"
          (fn [ax d o]
            (sv ye [1e10 -1e10])
            (defn update-ye [ye1]
              (sv (get ye 0) (min (first ye) (first ye1))
                  (get ye 1) (max (last ye) (last ye1))))
            (for [suffix (, "" "-src")]
              (for [(, ipg pg)
                    (enumerate (if (= suffix "")
                                   (, 2)
                                   (, 2
                                      (fn [np] np)
                                      (fn [np] (- np 2)))))]
                (when (and (not np-minus-2) (= ipg 2)) (continue))
                (assoc o :cdrglb "caas-node" :timeint "exact" :prefine 0
                       :measure :l2 :ref-ooa-2 False :ic (+ ic suffix)
                       :pg pg :pat-line (if (= suffix "")
                                            "-"
                                            (case/eq ipg [0  "--"] [1 "-."] [2 ":"])))
                (sv ye1 (p.plot ax d o))
                (update-ye ye1)
                (assoc o :timeint "interp" :prefine 5 :ref-ooa-2 (= suffix ""))
                (sv ye1 (p.plot ax d o))
                (update-ye ye1)))
            (pl.ylim (, (/ (first ye) 4) (nextpow10 (second ye)))))
          :ref-ooa-2 True)))

;;; filament diagnostic

(defn fig-filament [c d-all]
  (defn get-pat [nstepfac ne]
    (+ (get {10 "b" 20 "r." 40 "k"} ne)
       (get {1 "-" 5 "--"} nstepfac)))
  (sv p (AccFig) o (p.get-defaults c)
      nes (, 20 40)
      degs (nes->degstrs nes))
  (with [(pl-plot (, 10 4) (+ c.fig-dir "filament"))]
    (for [(, igrid grid) (enumerate (, :v :t))]
      (for [(, inp np) (enumerate c.nps)]
        (sv spi (inc (+ (* igrid (len c.nps)) inp))
            ax (pl.subplot 2 (len c.nps) spi))
        (for [(, ine ne) (enumerate nes)
              nstepfac (, 1 5)]
          (sv d (get d-all (if (= np 4) "exact" "interp") "nondivergent" nstepfac
                     "pcsl" (if (= np 4) "caas" "caas-node") "caas"
                     (if (= np 4) 0 5) np "cos" 1))
          (when (none? (geton d ne)) (continue))
          (sv thr (get d ne :L :thr)
              fil (get d ne :L (if (or (= grid :v) (= np 4)) :fil :me-fil)))
          (pl.plot thr fil (get-pat nstepfac ne)
                   :label (.format "{} {}" (get (:deg degs) ine)
                                   (nstepfac->word nstepfac)))
          (pl.xticks [0 0.2 0.4 0.6 0.8 1] :fontsize (:fs o))
          (pl.xlim (, 0.05 1))
          (pl.yticks (if (< np 8)
                         [0 20 40 60 80 100 120]
                         [50 60 70 80 90 100 110])
                     :fontsize (:fs o))
          (my-grid)
          (when (= spi (* 2 (len c.nps)))
            (pl.legend :loc "center left" :fontsize (:fs o) :frameon False))
          (pl.ylim (case/in np
                            [(, 9 12) (if (= grid :v) (, 75 115) (, 75 105))]
                            [(, 8) (, 50 115)]
                            [:else (, -5 125)]))
          (sv yl (pl.ylim)
              extra "")
          (when (> np 4) (sv extra (.format " {} grid"
                                            (if (= grid :v) "v" "tracer")))))
        (pl.text 0.1 (+ (first yl) (* 0.05 (npy.diff yl)))
                 (.format (if 0 "$n_p$ {}{}" "$n_p$ {} {}") np extra)
                 :fontsize (+ (:fs o) 2))))))

;;; mixing diagnostic

(defn triplot-read-dat [fname]
    (sv raw (with [f (open fname "rb")]
            (.read f))
      n (first (struct.unpack "i" (get raw (slice 0 4))))
      data (npy.zeros (* 2 n)))
  (for [i (range (* 2 n))]
    (sv os (+ 4 (* 4 i))
        (get data i) (first (struct.unpack "f" (cut raw os (+ os 4))))))
  (sv cb (cut data 0 n)
      ccb (cut data n (* 2 n)))
  (, cb ccb))

(defn triplot [cb ccb &optional [data None]]
  (defn triplot-curve [x]
    (+ 0.9 (* -0.8 x x)))
  (sv x (npy.linspace 0.1 1 100)
      lw 2)
  (pl.plot x (triplot-curve x) "k-" :lw lw)
  (sv tl (triplot-curve (last x))
      br (triplot-curve (first x)))
  (pl.plot x (* (triplot-curve (first x)) (npy.ones (len x))) "k-" :lw lw)
  (pl.plot (npy.ones (len x)) (npy.linspace tl br (len x)) "k-" :lw lw)
  (pl.plot x (npy.linspace br tl (len x)) "k-" :lw lw)
  (pl.plot cb ccb "r." :markersize 1)
  (pl.xlim 0.05 1.05)
  (pl.ylim 0.05 0.95)
  (sv t [0.2 0.4 0.6 0.8 1])
  (pl.xticks t []) (pl.yticks t [])
  (my-grid)
  ;(pl.axis "off")
  (defn text [x y txt sym me-sym]
    (pl.text x (+ y 0.1) txt)
    (sv dx 0.08)
    (pl.text (+ x dx) (+ y 0.1) (.format "{:1.2e} (v)" (get data sym)))
    (pl.text (+ x dx) y (if (in me-sym data)
                             (.format "{:1.2e}" (get data me-sym))
                             "")))
  (unless (none? data)
    (pl.text 0.1 0.91 (.format "$n_p$ {}, {} step" (get data 'np)
                               (nstepfac->word (get data 'nstepfac)))
             :fontsize 12)
    (sv x 0.1)
    (text x 0.32 "$l_r$" 'lr 'me-lr)
    (text x 0.11 "$l_u$" 'lu 'me-lu)))

(defn figs-mixing [c d]
  (sv p (AccFig) o (p.get-defaults c))
  (for [ne (, 20 40)]
    (with [(pl-plot (, 10 4) (+ c.fig-dir "mixing-ne" (str ne)) :format "png")]
      (for [(, instepfac nstepfac) (enumerate (, 1 5))
            (, inp np) (enumerate (:nps o))]
        (sv spi (inc (+ (* instepfac (len (:nps o))) inp))
            ax (pl.subplot 2 (len (:nps o)) spi)
            e (get d (if (= np 4) "exact" "interp") "nondivergent" nstepfac "pcsl"
                   (if (= np 4) "caas" "caas-node") "caas" (if (= np 4) 0 5) np
                   "cos" 1 ne :L)
            (, cb ccb) (triplot-read-dat (+ c.data-dir (:mixing-file e)))
            me-mixing (if (= np 4) :mixing :me-mixing)
            data {'np np 'nstepfac nstepfac
                  'lr (get e :mixing :lr) 'me-lr (get e me-mixing :lr)
                  'lu (get e :mixing :lu) 'me-lu (get e me-mixing :lu)})
        (unless (and (zero? (get e :mixing :lo)) (zero? (get e me-mixing :lo)))
          (prf "{}: lo {} {}" (:mixing-file e) (get e :mixing :lo)
               (get e me-mixing :lo)))
        (triplot cb ccb :data data)))))

;;; slotted cylinders images

(defn img-slo-filament [c d direc img-idx outname &optional nps nps-right]
  (svifn nps (, 4 4 6 8) nps-right nps)
  (sv degs (nes->degstrs (, 20 40))
      gap-for-colorbar (!= (len nps) (len nps-right))
      fs 7)
  (with [(pl-plot (, 7.3 (+ (if gap-for-colorbar 0 0.5) (len nps)))
                  outname :format "pdf" :tight False)]
    (sv spi 0 axs [])
    (for [(, inp np) (enumerate nps)
          (, ine ne) (enumerate (, 20 40))
          nstepfac (, 1 5)]
      (inc! spi)
      (when (and (= ne 40) (not (in np nps-right))) (continue))
      (sv c (% (dec spi) 4)
          r (// (dec spi) 4)
          w (/ 1 4)
          h (/ 1 (len nps))
          ax (pl.axes [(* w c) (- 1 (* h r)) (* 0.95 w) h]))
      (unless gap-for-colorbar (.append axs ax))
      (sv timeint (if (= np 4) "exact" "interp")
          cdrglb (if (and (= np 4) (= inp 0)) "caas" "caas-node")
          prefine (if (= np 4) 0 5)
          fname (+ direc "/" (.format "ne{}-np{}-nstep{}-{}-{}-pr{}.bin"
                                      ne np (* ne 6 nstepfac) timeint cdrglb prefine))
          img (try (get (read-slmmir-io-arrays fname) img-idx)
                   (except [e Exception] (print e "couldn't read" fname img-idx))))
      (unless (none? img)
        (print ne nstepfac np (/ (- 0.1 (npy.min img)) 0.1) (- (npy.max img) 1))
        (draw-slmmir-image img)
        (pl.text (- (last img.shape) 500) (- (first img.shape) 80)
                 (.format "$n_p$ {} {}" np (cdr-name cdrglb)) :fontsize fs)
        (defn write [x] (if (zero? x) "0" (.format "{:1.1e}" x)))
        (unless (none? d)
          (try
            (sv keys (, timeint "nondivergent" nstepfac "pcsl" cdrglb "caas" prefine
                        np "slo" 1 ne)
                e (get d #*(+ keys (, :Lerr))))
            (pl.text 10 15 (.format (+ "$l_2$ {:1.1e} $l_\infty$ {:1.1e}\n"
                                       "$\phi_{{min}}$ {:1.1e} $\phi_{{max}}$ {:1.1e}")
                                    (:l2 e) (:li e) (:phimin e) (:phimax e))
                     :fontsize fs)
            (except [] (print "no data for" timeint "nondivergent" nstepfac "pcsl"
                              cdrglb "caas" prefine np "slo" 1 ne))))
        (when (<= spi 4)
          (pl.text (/ (second img.shape) 2) (+ 30 (first img.shape))
                   (.format "{}, {} step" (get (:deg degs) ine)
                            (nstepfac->word nstepfac))
                   :ha "center" :fontsize (inc fs)))))
    (sv bdys (npy.linspace -0.05 1.05 23))
    (if gap-for-colorbar
        (do
          (sv c 2 r (+ 0.5 (- (len nps) 2))
              ax (pl.axes [(* w c) (- 1 (* h r)) (* 2 0.95 w) (* 0.2 h)])
              c (pl.colorbar :cax ax :orientation "horizontal" :aspect 30 :shrink 0.9
                             :ticks (npy.linspace 0 1.1 12)
                             :boundaries bdys)))
        (sv c (pl.colorbar :ax axs :orientation "horizontal" :aspect 50 :shrink 0.7
                           :ticks (npy.linspace 0 1.1 12) :pad 0.025
                           :boundaries bdys)))
    (c.ax.tick-params :labelsize (inc fs))))

(defn img-slocyl-slide [c d direc img-idxs outname &optional [ylabel False]]
  ;; 3 cols: IC above vertical colorbar, midpoint, endpoint
  (sv ne 20 deg (first (:deg (nes->degstrs (, 20)))) nstepfac 5 fs 7
      times (, "start" "middle" "end")
      w (/ 1 3)
      h (/ 1 (len c.nps)))
  (with [(pl-plot (, 5.5 (len c.nps)) outname :tight False)]
    (for [col (range 3)
          row (range (len c.nps))]
      (when (and ylabel (zero? col))
        ((. (pl.gcf) text)
         -0.01 (- 1 (* (- row 0.5) h))
         (if (zero? row) "Standard" "$p$-refined")
         :fontsize (+ fs 2) :rotation 90 :ha "right" :va "center"))
      (when (and (zero? col) (not (zero? row))) (continue))
      (sv ax (pl.axes [(* w col) (- 1 (* h row)) (* 0.95 w) h])
          np (nth c.nps row)
          timeint (if (= np 4) "exact" "interp")
          cdrglb (if (and (= np 4) (= inp 0)) "caas" "caas-node")
          prefine (if (= np 4) 0 5)
          fname (+ direc "/" (.format "ne{}-np{}-nstep{}-{}-{}-pr{}.bin"
                                      ne np (* ne 6 nstepfac) timeint cdrglb prefine))
          img (try (get (read-slmmir-io-arrays fname) (nth img-idxs col))
                   (except [e Exception]
                     (print e "couldn't read" fname (nth img-idxs col)))))
      (when (zero? col)
        (sv (get img (= img 1)) (- 1 1e-16)))
      (draw-slmmir-image img)
      (pl.text (- (last img.shape) 500) (- (first img.shape) 80)
               (.format "$n_p$ {} {}" np (cdr-name cdrglb)) :fontsize fs)
      (when (zero? row)
        (pl.text (/ (second img.shape) 2) (+ 30 (first img.shape))
                 (.format "{}, {} step, {}" deg (nstepfac->word nstepfac)
                          (nth times col))
                 :ha "center" :fontsize (inc fs))))
    (sv bdys (npy.linspace -0.05 1.05 23)
        col 0 row (dec (len c.nps)))
    (if (= (len c.nps) 2)
        (sv ax (pl.axes [(* w col) (- 1 (* 0.5 h row)) (* 0.95 w) (* 0.15 h)])
            cb (pl.colorbar :cax ax :orientation "horizontal" :aspect 30 :shrink 0.9
                            :ticks (npy.linspace 0 1 6)
                            :boundaries bdys))
        (sv col 0 row (dec (len c.nps))
            ax (pl.axes [(* w (+ col 0.35)) (- 1.02 (* h row)) (* 0.1 w)
                         (* (dec (len c.nps)) 0.95 h)])
            cb (pl.colorbar :cax ax :orientation "vertical" :aspect 30 :shrink 0.9
                            :ticks (npy.linspace 0 1 11)
                            :boundaries bdys)))
    (cb.ax.tick-params :labelsize (inc fs))))

;;; toy chemistry

(defn toychem-diagnostic-parse [fname]
  (sv txt (.split (readall fname) "\n")
      d {}
      skip 0 l2s [] lis [])
  (for [ln txt]
    (when (in "cmd> 1" ln) (break))
    (inc! skip))
  (for [ln (cut txt skip)]
    (cond [(in "cmd>" ln)
           (assert (or (empty? l2s) (= (len l2s) (* 576 10))))
           (sv cmd ln c (parse-cmd ln)
               key (+ (cmd->key-base c) (, (:ne c)))
               l2s [] lis [])
           (unless (none? (geton d #*key))
             ;; we collected pg<np data but only care about the pg=np data.
             (print "overwriting:\n" (get d #*(+ key (, :cmd)))))
           (assoc-nested d (+ key (, :cmd)) cmd)
           (assoc-nested d (+ key (, :l2)) l2s)
           (assoc-nested d (+ key (, :li)) lis)]
          [(in "C cycle" ln) (sv cyc (int (last (.split ln))))]
          [(= "C " (cut ln 0 2))]
          [(in "toy " ln)
           (sv pos (.find ln "toy")
               (, - masscons - l2 li) (sscanf (cut ln pos) "s,f,s,f,f"))
           (.append l2s l2)
           (.append lis li)]))
  d)

(defn fig-toychem-diagnostic [c d]
  (sv npa npy.array
      p (AccFig) o (p.get-defaults c))
  (with [(pl-plot (:figsize o) (+ c.fig-dir "toychem-diagnostic"))]
    (sv ax (pl.subplot 1 1 1))
    (for [npi (if 0 (, 8) (, 4.1 6 8 9 12))]
      (sv np (math.floor npi)
          e (get d (if (< npi 5) "exact" "interp") "nondivergent" 3 "pcsl"
                 (if (= npi 4) "caas" "caas-node") "caas" (if (< np 5) 0 5)
                 np np 30)
          l2 (npa (:l2 e))
          li (npa (:li e))
          n 576
          x (/ (npa (list (range (len l2)))) n)
          subset (list (range 0 (len l2) n))
          xsparse (/ (npa subset) n)
          clr (get c.npclrs np)
          mrk (get c.npmarks np))
      (pl.semilogy x l2 (+ clr "-") x li (+ clr ":")
                   xsparse (get l2 subset) (+ clr mrk)
                   xsparse (get li subset) (+ clr mrk)
                   :fillstyle "none" :lw (:lw o) :markersize (:markersize o)))
    (pl.ylabel "$\log_{10}$ Toy chemistry diagnostic" :fontsize (:fs o))
    (pl.xlabel "Cycles" :fontsize (:fs o))
    (svb (x (list (range 0 11)))
         (xt (cut x))
         #_((get xt 10) ""))
    (pl.xticks x xt :fontsize (:fs o))
    (sv y (npa (list (range -16 1))))
    (pl.yticks (** 10.0 y) y :fontsize (:fs o))
    (pl.title (+ "Toy chemistry diagnostic:\n"
                 "nondivergent flow, 1$^\circ$, 576 time steps/cycle\n"
                 "$p$-refinement, property preservation")
              :fontsize (:fs o))
    (p.legend ax (, (, "k-" "$l_2$") (, "k:" "$l_\infty$"))
              :o o :bbox (, 0 0.8))
    (pl.ylim (, (** 10.0 -16) (** 10.0 -3)))
    (pl.xlim (, -0.2 10.2))
    (my-grid)))

(defn img-toychem [c img-files outname &optional [diagnostic False]]
  (defn parse-filename [fname]
    (sv ode (cond [(in "nondiv" fname) "nondivergent flow"]
                  [(in "div" fname) "divergent flow"]
                  [:else "unknown-ode"])
        (, np pgtype pg) (cond [(in "np8pg8" fname) (, 8 "$n_f$" 8)]
                               [(in "np4pg0" fname) (, 4 "$n_p$" 4)])
        cdrtype (cond [(in "caas-node" fname) "CAAS-point"]
                      [(in "caas" fname) "CAAS-CAAS"])
        res (cond [(in "ne30" fname) "1$^\circ$"]
                  [:else "unknown-res"]))
    {:np np :pgtype pgtype :pg pg :dt "30min" :res res :ode ode :cdr cdrtype})
  (sv npa npy.array p (AccFig) o (p.get-defaults)
      ncol 2 nrow (// (+ (len img-files) ncol -1) ncol)
      vmin 0 vmax 4.5e-6 ticks (npy.linspace vmin vmax 10)
      ncolor (if diagnostic 8 18) fs 12
      xtbar 4e-6)
  (with [(pl-plot (, (* 3.65 ncol) (+ (* 2 nrow) 0.25)) outname :tight False)]
    (sv axs [])
    (for [spi (range (len img-files))]
      (sv c (% spi ncol) r (// spi ncol)
          w (/ 1 ncol) h (/ 1 nrow)
          ax (pl.axes [(* w c) (- 1 (* h r)) (* 0.95 w) h])
          img (read-slmmir-io-arrays (nth img-files spi))
          img (if diagnostic
                  (/ (- (+ (nth img (- (len img) 2))
                           (nth img (- (len img) 1)))
                        xtbar)
                     xtbar)
                  (nth img (- (len img) 2))))
      (sv imin (npy.min img) imax (npy.max img))
      ;; the toy chem paper draws images circle-shifted by 1/2 the image
      ;; relative to the bakeoff and test suite paper. do that here using
      ;; :switch-halves False if desired.
      (sv sh True)
      (if diagnostic
          (do
            (sv vmax (max (- imin) imax)
                thr (/ 0.3e-6 xtbar)
                vmax (if (> vmax thr) thr vmax)
                vmin (- vmax)
                ticks (npy.linspace vmin vmax 5))
            (draw-slmmir-image img :vmin vmin :vmax vmax :ncolor ncolor
                               :colorsym True :switch-halves sh))
          (draw-slmmir-image img :vmin vmin :vmax vmax :ncolor ncolor
                              :switch-halves sh))
      (sv d (parse-filename (nth img-files spi)))
      (pl.text 10 70
               (.format "$n_p$ {} {} {} {}"
                        (:np d) (:pgtype d) (:pg d) (:cdr d))
               :fontsize fs)
      (pl.text (- (last img.shape) 500) (- (first img.shape) 60)
               (.format "min {:8.1e} max {:8.1e}" imin imax)
               :fontsize (dec fs))
      (.append axs ax)
      (when diagnostic
        (sv c (pl.colorbar :ax ax :orientation "horizontal" :aspect 30 :shrink 0.8
                           :pad 0.05 :ticks ticks))
        (c.ax.tick-params :labelsize fs)
        (c.ax.set-xticklabels (lfor e ticks (.format "{:1.1e}"  e)))))
    (pl.text -20 (+ 15 (first img.shape))
             (.format "Toy chemistry: {}, $\Delta t$ {}, {}"
                      (:res d) (:dt d) (:ode d))
             :ha "center" :fontsize fs)
    (unless diagnostic
      (sv c (pl.colorbar :ax axs :orientation "horizontal" :aspect 50 :shrink 0.8
                         :pad 0.05 :ticks ticks))
      (c.ax.tick-params :labelsize fs)
      (c.ax.set-xticklabels (lfor e ticks (.format "{:1.1e}"  e))))))

;;; ISL comm footprint

(defn parse-footprint [fname]
  (sv d {}
      txt (.split (readall fname) "\n")
      c None cnt 0)
  (for [ln txt]
    (sv ln4 (cut ln 0 4))
    (cond [(= ln4 "cmd>")
           (dont unless (none? c)
             (print cmd)
             (print cnt))
           (sv cmd ln c (parse-cmd ln) cnt 0)]
          [(and (not (none? c)) (in "footprint>" ln))
           (sv pos (.find ln "t>")
               vals (sscanf (cut ln (+ pos 2)) "i,i,f,i"))
           (for [(, isym sym) (enumerate (, :min :median :mean :max))]
             (assoc-nested-append d (+ (cmd->key-base c) (, (:ne c) sym))
                                  (nth vals isym)))
           (inc! cnt)]))
  d)

(defn fig-comm-footprint [c d]
  (sv npa npy.array
      p (AccFig) o (p.get-defaults c))
  (assoc o :nps (, 4 8 12) :lw 1.5)
  (with [(pl-plot (:figsize o) (+ c.fig-dir "isl-footprint") :tight False)]
    (sv pat-line {:min ":" :median "-" :max "--"})
    (for [(, spi nstepfac) (enumerate (, 1 5))]
      (sv ax (pl.axes (, (/ spi 2) 0 0.41 0.8)))
      (for [np (:nps o)]
        (sv e (get d (if (= np 4) "exact" "interp") "nondivergent" nstepfac "pcsl"
                   (if (= np 4) "caas" "caas-node") "caas" (if (= np 4) 0 5)
                   np ; accidentally ran with pg; no effect but changes the dict nesting
                   np 30)
            x (* 12 (/ (npa (list (range (len (get e :min))))) (len (get e :min))))
            step (// (len x) 5)
            s (slice (// step 2) -1 step)
            xsparse (get x s))
        (for [sym (, :median :max)]
          (sv y (npa (get e sym))
              ysparse (get y s))
          (pl.plot [(first x) (last x)] (* (** np 2) (npa [1 1])) (+ (get c.npclrs np) ":"))
          (pl.plot x y
                   (+ (get c.npclrs np) (get pat-line sym))
                   xsparse ysparse
                   (+ (get c.npclrs np) (get c.npmarks np)) :fillstyle "none"))
        (my-grid)
        (pl.xticks (, 0 3 6 9 12) :fontsize (:fs o))
        (pl.yticks (npy.linspace 0 160 17) :fontsize (:fs o))
        (pl.ylim (, 0 160))
        (pl.xlabel "Days" :fontsize (:fs o))
        (pl.title (.format "{} time step" (get {1 "Long" 5 "Short"} nstepfac))
                  :fontsize (:fs o)))
      (cond [(= nstepfac 5)
             (p.legend ax (, (, (+ "k" (get pat-line :max)) "max")
                             (, (+ "k" (get pat-line :median)) "median")
                             (, (+ "k" ":") "$n_p^2$ reference"))
                       :o o :bbox (, 0 0.75) :nps-legend False)]
            [(= nstepfac 1)
             (p.legend ax (, ) :o o)])
      (when (= nstepfac 1)
        (pl.text 13 172
                 (+ "Islet: Number of transmitted scalars per tracer per element\n"
                    "Nondivergent flow")
                 :ha "center" :fontsize (:fs o))))))

;;; miscellaneous figs likely not to go in the paper

(defn img-slo-cyl-tracer-grid [c direc outname &optional [ne 10]]
  (defn make-fname [direc ne np nstepfac timeint cdrglb prefine]
    (+ direc "/" (.format "ne{}-np{}-nstep{}-{}-{}-pr{}.bin"
                          ne np (* ne 6 nstepfac) timeint cdrglb prefine)))
  (sv nstepfac 5)
  (sv degs (nes->degstrs (, 10))
      fs 7)
  (with [(pl-plot (, 8.25 4.5) outname)]
    (sv imgss [] nps [])
    (for [c (, (, 4 "exact" "caas" 0)
               (, 16 "interp" "caas-node" 5)
               (, 16 "interp" "caas-node" 1))]
      (.append nps (first c))
      (.append imgss (read-slmmir-io-arrays (make-fname direc ne (first c)
                                              nstepfac (second c) (nth c 2) (last c)))))
    (sv spi 0)
    (for [idx (, 1 3 5)
          (, i imgs) (enumerate imgss)]
      (sv np (nth nps i)
          img (nth imgs idx))
      (pl.subplot 3 3 (inc! spi))
      (when (= idx 1)
        (sv (get img (= img 1)) (- 1 1e-16)))
      (draw-slmmir-image img)
      (when (>= spi 7)
        (pl.text 10 25 (.format (+ "$n_e$ {} $n_p$ {} {} time step\n"
                                   "{}")
                                ne np (if (= nstepfac 5) "short" "long")
                                (case/in spi
                                         [(, 7 8) "on dynamics grid"]
                                         [:else "on tracer grid"]))
                 :fontsize fs)))))

;;; drivers

(when-inp ["acc-print-txt-table" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (when (in "stab-cmp" fname) (sv c.ics (, "gau")))
  (acc-print-txt-table c d))

(when-inp ["fig-stab-cmp" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (fig-stab-cmp c d))

(when-inp ["figs-acc" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (figs-acc c d))

(when-inp ["fig-midpoint" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (fig-acc-midpoint-check c d))

(when-inp ["fig-filament" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (fig-filament c d))

(when-inp ["figs-mixing" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (figs-mixing c d))

(when-inp ["img-filament" {:fname str :direc str}]
  (sv c (get-context)
      c.nps (, 4 6 8)
      d (acc-parse (+ c.data-dir fname)))
  (img-slo-filament c None (+ c.data-dir direc) 3 (+ c.fig-dir "slo-midpoint")
                    :nps (, 4 6 8 12) :nps-right (, 4 6 8))
  (img-slo-filament c d (+ c.data-dir direc) 5 (+ c.fig-dir "slo-finpoint")
                    :nps (, 4 6 8 12) :nps-right (, 4 6 8)))

(when-inp ["img-filament-slide" {:direc str}]
  (sv c (get-context)
      c.nps (, 4 6 8 12))
  (img-slocyl-slide c None (+ c.data-dir direc) (, 1 3 5)
                    (+ c.fig-dir "slo-imgs-slide"))
  (sv c.nps (, 4 12))
  (img-slocyl-slide c None (+ c.data-dir direc) (, 1 3 5)
                    (+ c.fig-dir "slo-imgs-slide-brief")
                    :ylabel True))

(when-inp ["fig-pg-mimic-src-term" {:fname str}]
  (sv c (get-context)
      d (acc-parse (+ c.data-dir fname)))
  (fig-acc-mimic-src-term-midpoint c d))

(when-inp ["fig-toychem-diagnostic" {:fname str}]
  (sv c (get-context)
      d (toychem-diagnostic-parse (+ c.data-dir fname)))
  (fig-toychem-diagnostic c d))

(when-inp ["fig-toychem-finpoint" {:direc str}]
  (sv c (get-context)
      d (+ c.data-dir direc "/")
      img-files [])
  (for [np (, 4 8)]
    (.append img-files (first (glob.glob (+ d (.format "*np{}*bin" np))))))
  (img-toychem c img-files (+ c.fig-dir "toychem-finpoint"))
  (img-toychem c img-files (+ c.fig-dir "toychem-finpoint-diagnostic")
               :diagnostic True))

(when-inp ["fig-comm-footprint" {:fname str}]
  (sv c (get-context)
      d (parse-footprint (+ c.data-dir fname)))
  (fig-comm-footprint c d))