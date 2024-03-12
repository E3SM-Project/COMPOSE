(require [amb3 [*]])
(import amb3 [amb3 [*]]
        [figsutils :as futils]
        [matplotlib :as mpl]
        [scipy.linalg :as linalg]
        math re poly islet)

(assoc matplotlib.rcParams "savefig.dpi" 300)
(do (pl-require-type1-fonts))

;;; tables

(defn read-last-array [ln]
  (sv arr (re.findall "{.*}" ln))
  (when (empty? arr) (return arr))
  (as-> (last arr) it
        (.replace it "," " ")
        (.replace it "{" "[")
        (.replace it "}" "]")
        (read-str it)
        (eval it)))

(defn parse-cpp-methods [fname]  
  (sv txt (.split (readall fname) "\n")
      d {})
  (for [ln txt]
    (cond [(in "xnodes: " ln)
           (sv xnodes (cond [(in "GLL" ln) 'gll]
                           [(in "Uniform" ln) 'uniform]))]
          [(in "case " ln)
           (sv np (int (first (re.findall "case (\d*):" ln))))]
          [(in "subnp[]" ln)
           (sv subnp (read-last-array ln)
               nodes [])]
          [(in "offst[]" ln)
           (sv offst (read-last-array ln)
               subtype 'offset)]
          [(and (not (in "nodes[]" ln))
                (not (empty? (do (sv arr (read-last-array ln))
                                 arr))))
           (sv subtype 'general)
           (.append nodes arr)]
          [(in "eval" ln)
           (assoc-nested d (, xnodes np)
                         (if (= subtype 'general)
                             {:subtype subtype :subnp subnp :nodes nodes}
                             {:subtype subtype :subnp subnp :offst offst}))]))
  d)

(defn make-array-str [a] (+ "\{" (str.join ", " (lfor e a (str e))) "\}"))

(defn make-subsets-str [ns]
  (sv s "\n\\begin{tabular}{l}\n\!\!\!nodal subsets \\\\\n\{")
  (for-first-last
    [n ns]
    (+= s
        (if first? "" "\\phantom{\{}")
        (make-array-str n)
        (if last? "\}\n" ",\\\\\n")))
  (+ s "\end{tabular}"))

(defn write-methods-latex-tables [d xnodes fname]
  (sv col-hdrs (, "$\\numg$" "OOA" "$\\npsub$" "Supports")
      d (get d xnodes)
      nps (list (.keys d)))
  (with [f (open fname "w")]
    (defn write [s] (.write f s))
    (write "\\begin{center}\n")
    (write "\\begin{tabular}{rcll}\n")
    (write (+ (str.join " & " col-hdrs) " \\\\\n\\hline\n"))
    (when (= xnodes 'gll)
      (write "4 & 2 & see text & see text \\\\\n\\hline\n"))
    (for [np nps]
      (sv e (get d np))
      (svb (np-str (str np))
           (nodes-str (case/eq xnodes ['gll "GLL"] ['uniform "Uniform"]))
           (ooa-str (str (dec (min (:subnp e)))))
           (npsub-str (make-array-str (:subnp e)))
           (subsets-str (if (= (:subtype e) 'offset)
                            (+ "offsets " (make-array-str (:offst e)))
                            (make-subsets-str (:nodes e)))))
      (write (+ (str.join " & " (, np-str ooa-str npsub-str subsets-str))
                " \\\\\n"))
      (unless (= np (last nps)) (write "\\hline\n")))
    (write "\\end{tabular}\n")
    (write "\\end{center}\n")))

;;; illustrations

(defn draw-gll-dots [np marker &optional [dx 0] [dy 0] [markersize 12]
                     fillstyle [alpha 1]]
  (sv xgll (col-vec (poly.get-gll-x np))
      e (npy.ones (, np 1))
      X (npy.dot e (.transpose xgll))
      Y (npy.dot xgll (.transpose e)))
  (pl.plot (+ dx X) (+ dy Y) marker :markersize markersize :fillstyle fillstyle))

(sv element-blue "#3838FF")

(defn illustrate-grids [img-fname]
  (sv nf 6 elc "k")
  (with [(pl-plot (, 5 5) img-fname :format "pdf")]
    (pl.plot [-1 1] [1 1] elc [-1 1] [-1 -1] elc [1 1] [-1 1] elc [-1 -1] [-1 1] elc
             :linewidth 2 :color element-blue)
    (sv lw 2
        line "g--"
        d 0.99 nd (- d))
    (for [it (range 2)]
      (for [i (range (inc nf))]
        (sv x (as-> (- (* 2 (/ i nf)) 1) x
                    (if (= x -1) nd
                        (= x 1) d
                        x)))
        (if (zero? it)
            (pl.plot [x x] [nd d] line :linewidth lw)
            (pl.plot [nd d] [x x] line :linewidth lw))))
    (draw-gll-dots 4 (+ elc "o"))
    (draw-gll-dots 8 "ro" :markersize 7)
    (pl.axis "equal")
    (pl.axis "off")))

(defn draw-np4-schematic [img-fname]
  (sv np 4 fs 12
      x-gll (poly.get-gll-x np)
      nx 256
      x (npy.linspace -1 1 nx)
      isl (islet.Islet)
      (, yn yo yb) (lfor m (, 0 2 1) (.transpose (isl.eval m np x)))
      clrs "kgrb"
      d 1.04)
  (with [(pl-plot (, 4 4) img-fname :tight True)]
    (pl.axes [0 0.6 1 0.4])
    (sv xsub (npy.linspace (nth x-gll 2) 1 (// nx 4))
        xscale (- (* 2 (/ (- xsub (nth x-gll 2)) (- 1 (nth x-gll 2)))) 1)
        alpha (poly.eval-lagrange-poly (poly.get-gll-x 3) (, 0 0.306 1) xscale))
    (pl.plot (- (get xsub (s-all-rev))) (get alpha (s-all-rev)) "k-."
             xsub alpha "k-.")
    (do
      (pl.plot -2 -2 "k:"  :label "Natural")
      (pl.plot -2 -2 "k--" :label "Offset nodal subset: $\{3,4\}$, $\{0,0\}$")
      (pl.plot -2 -2 "k-"  :label "Optimized")
      (pl.plot -2 -2 "k-." :label "Convex combination parameter")
      (pl.legend :loc "center" :fontsize fs :bbox-to-anchor (, 0.49 0.64)
                 :frameon False))
    (my-grid)
    (axis-tight-pad :pad 0)
    (pl.xlim (, (- d) d))
    (pl.ylim (, -0.05 1.05))
    (pl.yticks (, 0 0.5 1) :fontsize fs)
    (pl.plot (nth x-gll 2) 0 "ko" 1 1 "ko")
    (pl.text 0.57 1 "Fully $n_p$ 4" :fontsize fs :va "top")
    (pl.text 0.05 0 "Fully $n_p$ 3" :fontsize fs)
    (pl.title "Optimized $n_p$ 4 basis" :fontsize fs)
    (sv xt (npy.linspace -1 1 9))
    (pl.xticks xt (lfor e (range (len xt)) "") :fontsize fs)
    (pl.axes [0 0 1 0.57])
    (for [i (range (dec np) -1 -1)]
      (sv c (nth clrs i))
      (pl.plot x (get yn i) (+ c ":"))
      (pl.plot x (get yo i) (+ c "--"))
      (pl.plot x (get yb i) (+ c "-")))
    (my-grid)
    (axis-tight-pad :pad 0)
    (pl.xlim (, (- d) d))
    (pl.xticks xt :fontsize fs)
    (pl.yticks (npy.linspace -0.2 1 7) :fontsize fs)
    (pl.xlabel "Reference-element coordinate" :fontsize fs)))

(defn draw-islet-timestep-schematic [img-fname &optional [npv 4] [npt 6]]
  (defn make-gll-grid [np]
    (sv x (poly.get-gll-x np)
        ones (npy.ones np)
        X (npy.dot (col-vec x) (row-vec ones))
        Y (npy.dot (col-vec ones) (row-vec x)))
    (, X Y))

  (defn rotate [X Y theta]
    (sv c (npy.cos theta) s (npy.sin theta))
    (, (+ (* (- s) Y) (* c X))
       (+ (* c Y) (* s X))))
  (defn scale [X Y p]
    (sv rad (npy.sqrt (+ (** X 2) (** Y 2)))
        s (+ 1 (* p (- rad 1))))
    (, (* s X) (* s Y)))
  (defn perturb [X Y p]
    (sv np (first X.shape)
        r (fn [] (as-> (npy.random.rand np np) x
                       (* 2 x)
                       (- x 1)
                       (* p x))))
    (, (+ X (r)) (+ Y (r))))
  (defn get-v-pts [np]
    (sv (, X0 Y0) (make-gll-grid np))
    (+= X0 0.3)
    (+= Y0 -0.2)
    (sv (, X0 Y0) (as-> (, X0 Y0) xy
                        (scale #*xy 0.3)
                        (perturb #*xy 0.05) xy
                        (rotate #*xy 0.2)))
    (+= X0 0.95)
    (+= Y0 -0.1)
    (, X0 Y0))

  (defn draw-v-elem [ax X Y npt v-marker v-marker-sz t-marker t-marker-sz]
    (sv vgll (poly.get-gll-x npv)
        f (fn [xi idxs]
            (sv x (get X idxs) y (get Y idxs))
            (-> (zip (poly.eval-lagrange-poly vgll x xi)
                     (poly.eval-lagrange-poly vgll y xi))
                (list)))
        g (fn [x]
            (+ (f x (, (s-all) -1))
               (f x (, -1 (s-all-rev)))
               (f x (, (s-all-rev) 0))
               (f x (, 0 (s-all)))))
        x (npy.linspace -1 1 51)
        vtxs (g x))
    (sv p (mpl.patches.Polygon vtxs :facecolor "w" :fill True :alpha 0.7
                               :zorder 10))
    (.add-patch ax p)
    (sv vtxs (npy.array (unzip vtxs)))
    (pl.plot (get vtxs 0) (get vtxs 1) "--" :color element-blue :zorder 11 :lw 3)
    (pl.plot X Y v-marker :markersize 14 :zorder 12
             :markerfacecolor "w" :markeredgecolor "k")
    (sv tgll (poly.get-gll-x npt)
        t-pts [])
    (for [i (range npt) j (range npt)]
      (sv vi (poly.eval-lagrange-poly-basis vgll (nth tgll i))
          vj (poly.eval-lagrange-poly-basis vgll (nth tgll j))
          basis (npy.dot (col-vec vi) (row-vec vj)))
      (.append t-pts (, (npy.sum (* basis X)) (npy.sum (* basis Y)))))
    (sv t-pts (npy.array (unzip t-pts))
        t-in-mask (npy.logical-and (>= (get t-pts 0) 1)
                                   (<= (get t-pts 1) 0))
        t-out-mask (npy.logical-not t-in-mask))
    (for [(, mask marker sz) (, (, t-in-mask "g*" 12)
                                (, t-out-mask "go" 8))]
      (pl.plot (get t-pts 0 mask) (get t-pts 1 mask) marker
               :markersize sz :zorder 13)))

  (defn plot-basis [np]
    (sv isl (islet.Islet)
        gll-best 1
        pats {0 "--" gll-best "-"}
        clrs "krbcgm"
        x (npy.linspace -1 1 213)
        fs 18)
    (for [method (, 0 gll-best)]
      (sv y (.transpose (isl.eval method np x))
          pat (get pats method))
      (for [i (range np)]
        (pl.plot x (get y i) (+ (nth clrs (% i (len clrs))) pat)
                 :label (if (zero? i)
                            (get {0 "Natural GLL" gll-best "Islet GLL"}
                                 method))))
      (my-grid)
      (sv d 1.04)
      (pl.yticks (, 0 1) :fontsize fs)
      (pl.xticks (, -1 0 1) :fontsize fs)
      (pl.xlim (, (- d) d))
      (pl.ylim (, -0.25 1.05)))
    (pl.xlabel "Reference-element coordinate" :fontsize fs)
    (pl.ylabel "Basis function value" :fontsize fs)
    (pl.text -1.03 1.1 (.format "Basis functions, $n_p$ = {}" np) :fontsize fs)
    (pl.figlegend :loc (, 0.42 0.36) :fontsize (- fs 1) :ncol 2))

  (sv nelx 3 nely 2
      f (fn [n] (- (/ n 2) 0.5))
      hnelx (f nelx) hnely (f nely)
      elw 2 f (fn [n] (* elw (- (/ n 2))))
      x0 (f nelx) y0 (f nely)
      (, X0 Y0) (get-v-pts npv)
      one (npy.ones 2))

  (with [(pl-plot (, 8 8) img-fname)]
    (sv ax (pl.axes (, 0 0.35 1 0.65))
        lw 3)
    (for [i (range (inc nelx))]
      (pl.plot (* (+ x0 (* i elw)) one) [y0 (+ y0 (* elw nely))] "-"
               :color element-blue :lw lw))
    (for [i (range (inc nely))]
      (pl.plot [x0 (+ x0 (* elw nelx))] (* (+ y0 (* i elw)) one) "-"
               :color element-blue :lw lw))
    (for [i (range nelx) j (range (dec nely) -1 -1)]
      (sv pats (if (and (= i 2) (= j 0))
                   (, (, "ko" "full" 11) (, "ro" "full" 6))
                   (, (, "ro" "full" 8))))
      (for [pat pats]
        (draw-gll-dots npt (first pat)
                       :dx (* elw (- i hnelx)) :dy (* elw (- j hnely))
                       :markersize (last pat) :fillstyle (second pat)
                       :alpha 0.5)))
    (draw-gll-dots npv "ko" :markersize 12 :dx -2 :dy 1)
    (draw-v-elem ax X0 Y0 npt "ko" 12 "ro" 8)
    (pl.axis "tight")
    (pl.axis "equal")
    (pl.axis "off")
    (pl.text -3.045 -2.3 "(a)" :fontsize 18)
    (sv x0 0.9 x1 3.1 x1m 2.9 y0 -2.1 y1 -1.9 y2 0.1)
    (pl.plot [x0 x1 x1 x0 x0] [y0 y0 y1 y1 y0] "r--")
    (pl.plot [x1m x1 x1 x1m x1m] [y0 y0 y2 y2 y0] "g--")
    (pl.axes (, 0.1 0 0.85 0.3))
    (plot-basis 6)
    (pl.text -1.19 -0.5 "(b)" :fontsize 18)))

;;; utils for search

(defn parse-search-list [fname]
  (sv txt (.split (readall fname) "\n")
      c [])
  (for [ln txt]
    (when (or (< (len ln) 100) (!= (cut ln 0 5) "meam1")) (continue))
    (sv b (futils.parse-search-basis-line ln))
    (.append c b))
  c)

(defn uniquify-search-list [c]
  (sv uniq (set) c-uniq [])
  (for [b c]
    (sv s (futils.Nodes->string (:np b) (:nodes b)))
    (when (in s uniq) (continue))
    (.add uniq s)
    (.append c-uniq b))
  c-uniq)

(defn get-slmmir-builtin [np]
  (get
    {4  None
     5  None
     6  [[0 1 2 3 4] [ 0 1 2 3 5] [0 1 2 3 4 5]]
     7  (futils.offst->Nodes (, 5 5 6) (, 0 0 0))
     8  (futils.offst->Nodes (, 6 6 7 6) (, 0 0 0 1))
     9  [[0 1 2 3 4 5 8] [0 1 2 3 4 5 7 8] [0 1 2 3 4 5 6 8] [1 2 3 4 5 6 7]]
     10 (futils.offst->Nodes (, 7 7 7 8 8) (, 0 0 0 0 1))
     11 (futils.offst->Nodes (, 8 9 8 9 8) (, 0 0 0 0 1))
     12 (futils.offst->Nodes (, 9 9 10 10 9 10) (, 0 0 0 0 1 1))
     13 (futils.offst->Nodes (, 10 10 10 10 11 10) (, 0 0 0 0 0 1))}
    np))

(defn nodes= [a-np a-nodes b-np b-nodes]
  (unless (and (= a-np b-np)
               (= (len a-nodes) (len b-nodes)))
    (return False))
  (for [i (range (len a-nodes))]
    (sv npa npy.array ai (npa (nth a-nodes i)) bi (npa (nth b-nodes i)))
    (unless (and (= (len ai) (len bi)) (npy.all (= ai bi)))
      (return False)))
  True)

;;; run slmmir on a bunch of formulas

(defn write-slmmir-script [blns script-fname]
  (svb (cmdstr
         (fn [basis]
           (+ "KMP_AFFINITY=balanced OMP_NUM_THREADS=48 $exe "
              "-method pcsl -ode {ode:s} -ic gaussianhills -ic cosinebells -ic slottedcylinders "
              "-we 0 -rit -dmc eh -mono {mono:s} -lim {lim:s} -nsteps {nstep:d} -T 12 -ne {ne:d} "
              "-np {np:d} -timeint {timeint:s} -prefine {prefine:d} -d2c "
              (if basis "-basis \"{basis:s}\" " "")
              "|& grep \"^C \"")))
       (ode "nondivergent") (ne 20) (nstep (* 20 6))
       (lims (, (, "caas-node" "caas") (, "none" "none"))))
  (with [f (open script-fname "w")]
    (f.write "exe=\n")
    (sv ctr 0)
    (defn write1 [b lim str-basis]
      (sv cmd (.format (cmdstr str-basis) :ode ode :mono (first lim)
                       :lim (last lim) :nstep nstep :ne ne :np (:np b)
                       :timeint (if (= (:np b) 4) "exact" "interp")
                       :prefine (if (= (:np b) 4) 0 5)
                       :basis (if str-basis (futils.Nodes->string (:np b) (:nodes b)))))
      (f.write (.format "echo 'line> {}'\n" (if str-basis
                                                (:txt b)
                                                (.format "builtin np {}" (:np b)))))
      (f.write (.format "echo 'cmd> {} {}'\n" ctr cmd))
      (f.write (.format "{}\n" cmd)))
    (dont for [np (range 4 13) lim lims]
      (write1 {:np np} lim False)
      (inc! ctr))
    (for [b blns lim lims]
      (write1 b lim True)
      (inc! ctr))))

(defn parse-slmmir-output [fname &optional d [lebesgue False]]
  (import islet)
  (svifn d {})
  (sv txt (.split (readall fname) "\n")
      isl (islet.Islet))
  (for [ln txt]
    (cond [(in "line>" ln)
           (sv b (futils.parse-search-basis-line (cut ln 6)))
           (when lebesgue
             (sv bstr (futils.Nodes->string (:np b) (:nodes b))
                 npm (isl.calc-xnodes-metrics-from-basis-string bstr)
                 npa npy.array)
             (assert (< (reldif (npa (:npm b)) (npa npm)) 1e-2))
             (assoc b :lebesgue (isl.calc-lebesgue-consts-from-basis-string bstr))
             (sv strp (fn [t]
                        (sv s "")
                        (for [e t] (+= s (.format " {:1.2e}" e)))
                        s))
             (dont prf (+ (strp (:npm b)) " |" (strp (:lebesgue b)))))]
          [(in "cmd>" ln)
           (sv cmd ln c (futils.parse-cmd ln)
               cls {})]
          [(in "C cycle" ln)
           (sv cyc (int (last (.split ln))))]
          [(or (in "C PASS" ln) (in "C FAIL" ln))
           (assoc-nested-append d (+ (futils.cmd->key-base c)
                                     (, cyc (:ne c)))
                                {:b b :cls cls})]
          [(and (> (len ln) 10) (= (cut ln 0 2) "C "))
           (sv cl (futils.parse-C ln))
           (assoc cls (:ic cl) cl)]))
  d)

(defn basic-key [np &optional ode nstepfac prop-preserve cyc ne]
  (svifn ode "nondivergent" nstepfac 1 cyc 1 ne 20)
  (, (if (= np 4) "exact" "interp") ode nstepfac "pcsl"
     (if prop-preserve "caas-node" "none") (if prop-preserve "caas" "none")
     (if (= np 4) 0 5) np cyc ne))

(defn get-tick-labels [x]
  (sv xt [])
  (for [e x]
    (sv v (npy.log10 e))
    (.append xt (if (= v (int v)) (int v) (.format "{:1.1f}" v))))
  xt)

(defn set-semilogy-ticks [&optional fs]
  (sv y (first (pl.yticks)))
  (pl.yticks y (get-tick-labels y) :fontsize fs))

(defn set-loglog-ticks [&optional fs]
  (sv x (first (pl.xticks)))
  (pl.xticks x (get-tick-labels x) :fontsize fs)
  (sv y (first (pl.yticks)))
  (pl.yticks y (get-tick-labels y) :fontsize fs))

(defn norm->str [n]
  (get {:l1 "$l_1$" :l2 "$l_2$" :li "$l_{\infty}$"} n))

(defn plot1-slmmir-vs-heuristic [c d nps prop-preserve ic norm pum-thr lebesgue
                                 &optional [jcp False]]
  (sv npa npy.array fs 11
      plot (if lebesgue pl.semilogy pl.loglog))
  (for [np nps]
    (sv es (get d #*(basic-key np :prop-preserve prop-preserve)))
    (when (none? es) (continue))
    (sv x [] y [])
    (for [e es]
      (sv cls (:cls e) b (:b e))
      (when (> (:pum b) pum-thr) (continue))
      (.append x (nth (get b (if lebesgue :lebesgue :npm))
                      (get {:l1 0 :l2 1 :li 2} norm)))
      (.append y (get cls ic norm))
      (when (nodes= (:np b) (:nodes b) (:np b) (get-slmmir-builtin (:np b)))
        (plot (last x) (last y) "ro" :fillstyle "none" :markersize 12
              :zorder 20)))
    (plot (npa x) (npa y) (+ (get c.npclrs np) (get c.npmarks np))
          :fillstyle "none" :zorder (- 20 np) :label (.format "$n_p$ {}" np)))
  (if lebesgue
    (do (set-semilogy-ticks :fs fs)
        (pl.xlabel (+ (norm->str norm) " heuristic") :fontsize fs))
    (do (set-loglog-ticks :fs fs)
        (pl.xlabel (+ "$\log_{10}$ " (norm->str norm) " heuristic")
                   :fontsize fs)))
  (my-grid)
  (pl.title (.format (+ "Test problem vs.~heuristic:\n"
                        "nondivergent flow, "
                        (if jcp "$n_e = 20$, " "1.5$^\circ$, ")
                        "{}, "
                        (if jcp "120 steps\n" "long steps\n")
                        "$p$-refinement, {}")
                     (futils.ic-short2long ic)
                     (+ (if prop-preserve "" "no ") "property preservation"))
            :fontsize fs)
  (pl.legend :loc "best" :fontsize fs)
  (pl.axis "tight")
  (pl.ylabel (+ "$\log_{10}$ " (norm->str norm) " relative error")
             :fontsize fs))

(defn plot-slmmir-vs-heuristic [c d img-fname
                                &optional nps prop-preserve ic norm pum-thr
                                lebesgue]
  (svifn nps (list (range 4 14)) prop-preserve False ic "gau" norm :l2
         pum-thr 1e-6 lebesgue False)
  (sv npa npy.array fs 11
      plot (if lebesgue pl.semilogy pl.loglog))
  (print img-fname)
  (with [(pl-plot (, 4 4.2) img-fname :format "pdf")]
    (plot1-slmmir-vs-heuristic c d nps prop-preserve ic norm pum-thr lebesgue)))

(defn plot-slmmir-vs-heuristic-ab [c d img-fname norm-pp-ic-tuples
                                   &optional nps pum-thr]
  (svifn nps (list (range 4 14)) pum-thr 1e-6)
  (sv lebesgue False)
  (sv npa npy.array fs 11
      plot pl.loglog)
  (with [(pl-plot (, 7.8 4.2) img-fname :format "pdf")]
    (for [i (range 2)]
      (sv t (nth norm-pp-ic-tuples i))
      (pl.subplot 1 2 (inc i))
      (plot1-slmmir-vs-heuristic c d nps (nth t 1) (nth t 2) (nth t 0) pum-thr
                                 lebesgue :jcp True)
      (sv f (pl.gcf))
      (f.text (nth (, 0.02 0.52) i) 0.05 (nth (, "(a)" "(b)") i)
              :fontsize fs))))

;;; pum vs perturb

(defn parse-pum-vs-perturb [fname]
  (sv txt (.split (readall fname) "\n")
      perturb None
      d {})
  (for [ln txt]
    (cond [(= ">>> " (cut ln 0 4))
           (sv (, - basis np) (sscanf ln "s,s,i"))]
          [(= ">> " (cut ln 0 3))
           (sv pos (.find ln "[")
               arr (eval (read-str (cut ln pos))))
           (if (none? perturb)
               (sv perturb arr)
               (do (sv meam1 arr)
                   (assoc-nested d (, basis np) (, perturb meam1))))]))
  d)

(defn pum-make-reference-slope-triangle
  [x-span y-span slope pattern
   &optional [opposite False] [kwargs-plot None] [kwargs-text None]]
    (assert (= 2 (len x-span)))
    (assert (= 2 (len y-span)))
    (svifn kwargs-plot {})
    (svifn kwargs-text {})
    (sv (, x0 x1) x-span (, y0 y1) y-span
        dx (- x1 x0) dy (- y1 y0)
        x [x0 x1 x1 x0]
        y [y1 y1 y0 y1])
    (pl.plot #*[x y pattern] #**kwargs-plot)
    (pl.text #*[(* 0.7 x1)
                (* 0.05 y0)
                (str slope)]
             #**kwargs-text))

(defn plot-pum-vs-perturb [c d fname]
  (sv f identity fs 11)
  (with [(pl-plot (, 4 4) fname)]
    (for [basis (, "gll_best")]
      (sv uni (in "uni" basis))
      (for [np (range (if uni 8 6) 14)]
        (sv e (geton d basis np))
        (when (none? e) (continue))
        (pl.loglog (f (first e)) (f (second e))
                   (+ (get c.npclrs np) (if uni "--" "-"))
                   :label (+ "$n_p$ " (str np) (if uni "'" "")))))
    (sv f 1.7 slope 4)
    (pum-make-reference-slope-triangle
     [(/ 2.8e-3 f) (* 2.8e-2 f)]
     [(* 3e-8 (** f slope)) (/ 3e-12 (** f slope))]
     slope "k-"
     :kwargs-text {"fontsize" (+ fs 4)})
    (set-loglog-ticks)
    (my-grid)
    (pl.legend :loc "upper left" :fontsize (dec fs) :ncol 1)
    (pl.xlabel "$\log_{10}$ Element size relative random perturbation $\delta$"
               :fontsize fs)
    (pl.ylabel "$\log_{10}$ (max $|\lambda|$ - 1)" :fontsize fs)
    (pl.title "Perturbed uniform grid metric" :fontsize fs)
    (pl.xlim (, 1e-4 0.1))
    (pl.ylim (, 1e-14 0.1))))

;;; meam1 and pum vs dx

(defn parse-meam1-sweep [fname]
  (sv txt (.split (readall fname) "\n")
      methods (, "gll_natural" "gll_best" "uniform_offset_nodal_subset")
      fst True
      d {})
  (defn fill []
    (print method (len dxs) (len meam1s))
    (assoc d method (, dxs meam1s np)))
  (for [ln txt]
    (sv toks (.split ln))
    (cond [(in (first toks) methods)
           (unless fst (fill))
           (sv method (first toks)
               np (int (second toks))
               fst False
               dxs [] meam1s [])]
          [(= (len toks) 2)
           (sv (, dx meam1) (sscanf ln "f,f"))
           (.append dxs dx)
           (.append meam1s meam1)]))
  (fill)
  d)

(defn parse-pum-sweep [method-fnames]
  (sv d {})
  (for [(, method fname) method-fnames]
    (sv txt (.split (readall fname) "\n")
        dxs [] meam1s [] skip 1)
    (for [ln txt]
      (when (in "final" ln) (break))
      (inc! skip))
    (for [ln (cut txt skip)]
      (sv toks (.split ln))
      (unless (= (len toks) 4) (break))
      (.append dxs (float (nth toks 1)))
      (.append meam1s (float (nth toks 3))))
    (assoc d method (, dxs meam1s)))
  d)

(defn plot-meam1-and-pum-vs-dx [c dmeam1 dpum fname]
  (defn symx [x]
    (sv x (npy.array x))
    (npy.concatenate (, x (- 1 (get x (s-all-rev))))))
  (defn symy [y]
    (sv y (npy.array y))
    (npy.concatenate (, y (get y (s-all-rev)))))
  (sv fs 11
      ms (, "gll_natural" "uniform_offset_nodal_subset" "gll_best")
      clrs {(nth ms 0) "r" (nth ms 1) "g" (nth ms 2) "k"}
      mrks {(nth ms 0) "." (nth ms 1) "x" (nth ms 2) "o"}
      lbls {(nth ms 0) "GLL natural" (nth ms 1) "Uniform offset nodal subset"
            (nth ms 2) "GLL offset nodal subset"})
  (with [(pl-plot (, 4 4) fname)]
    (for [method ms]
      (sv e (get dpum method))
      (pl.semilogy (first e) (second e)
                   (+ (get clrs method) (get mrks method))
                   :label (get lbls method))
      (sv e (get dmeam1 method))
      (pl.semilogy (symx (first e)) (symy (second e))
                   (+ (get clrs method) "-")))
    (sv np (last e))
    (for [i (range (// np 2) np)]
      (sv p (/ i (dec np)))
      (pl.semilogy [p p] [1e-15 1] (+ (get clrs (nth ms 1)) ":") :zorder -1))
    (set-semilogy-ticks)
    (pl.ylim (, 1e-16 1))
    (pl.ylabel "$\log_{10}$ (max $|\lambda|$ - 1)" :fontsize fs)
    (pl.xlabel "Translation distance $\Delta x$ relative to element size 1" :fontsize fs)
    (pl.title (.format "$n_p$ {} methods" np) :fontsize fs)
    (pl.legend :loc (, 0.02 0.15) :fontsize fs)
    (pl.xticks (npy.linspace 0 1 11))
    (my-grid)))

;;; drivers

(when-inp ["dev-parse"]
  (sv fname "data/mar21/search-0.txt"
      blns (parse-search-list fname))
  (print (len blns))
  (for [b blns]
    (when (and (= (:np b) 13) (= (:type b) :offst-nodal-subset) (< (:pum b) 1e-6))
      (print b) (break)))
  (for [b blns]
    (when (and (= (:np b) 8) (= (:type b) :nodal-subset) (< (:pum b) 1e-7))
      (print b) (break))))

(when-inp ["write-slmmir-script"]
  (sv script-fname "../../slmm/meas/run-slmmir-on-basis-lines.sh"
      data-fnames (, "search-0.txt"
                     "search-findnodal_given_bestosn-0.txt"
                     "search-findnodal_given_bestosn-1.txt"
                     "search-findnodal_given_bestosn-2.txt"
                     "search-findnodal_given_bestosn-3.txt"
                     "search-findnodal_given_bestosn-7.txt")
      blns [])
  (for [fname data-fnames]
    (.extend blns (parse-search-list (+ "data/mar21/" fname))))
  (sv blns (uniquify-search-list blns))
  (write-slmmir-script blns script-fname))

(when-inp ["plot-slmmir-vs-heuristic"]
  (sv fnames (, "slmmir-on-basis-lines-2.txt")
      c (futils.get-context)
      lebesgue False
      d {})
  (for [fname fnames]
    (sv d (parse-slmmir-output (+ "data/mar21/" fname) :d d :lebesgue lebesgue)))
  (for [(, norm pp ic) (, (, :l2 False "gau") (, :l2 True "cos"))]
    (plot-slmmir-vs-heuristic
      c d (.format (+ "{}slmmir-vs-heuristic-{}-{}-{}" (if lebesgue "-leb" ""))
                   c.fig-dir ic (if pp "pp" "nopp") (name norm))
      :nps [6 7 8 9 10]
      :prop-preserve pp :norm norm :ic ic :lebesgue lebesgue)))

(when-inp ["plot-slmmir-vs-heuristic-ab"]
  (sv fnames (, "slmmir-on-basis-lines-2.txt")
      c (futils.get-context)
      lebesgue False
      d {})
  (for [fname fnames]
    (sv d (parse-slmmir-output (+ "data/mar21/" fname) :d d :lebesgue lebesgue)))
  (plot-slmmir-vs-heuristic-ab
   c d (+ c.fig-dir "slmmir-vs-heuristic-ab")
   (, (, :l2 False "gau") (, :l2 True "cos"))
   :nps [6 7 8 9 10]))

(when-inp ["tables"]
  (for [xnodes (, 'gll)]
    (write-methods-latex-tables (parse-cpp-methods "islet-methods.txt")
                                xnodes (.format "figs/methods-table-{}.tex"
                                                (name xnodes)))))

(when-inp ["pum-vs-perturb"]
  (sv fname "pum_perturb_plot-041021.txt"
      d (parse-pum-vs-perturb (+ "data/mar21/" fname))
      c (futils.get-context))
  (plot-pum-vs-perturb c d (+ c.fig-dir "pum-vs-perturb")))

(when-inp ["meam1-and-pum-vs-dx"]
  (sv c (futils.get-context)
      data-dir "data/mar21/"
      meam1-fname (+ data-dir "run_meam1_sweep-np8.txt")
      method-fnames (zip (, "gll_best" "gll_natural" "uniform_offset_nodal_subset")
                         (lfor fname (, "pum_sweep-np8-gll_best.txt"
                                        "pum_sweep-np8-gll_natural.txt"
                                        "pum_sweep-np8-uni.txt")
                               (+ data-dir fname)))
      d (parse-meam1-sweep meam1-fname)
      d-pum (parse-pum-sweep method-fnames))
  (plot-meam1-and-pum-vs-dx c d d-pum (+ c.fig-dir "meam1-and-pum-vs-dx")))

(when-inp ["np4-schematic"]
  (sv c (futils.get-context))
  (draw-np4-schematic (+ c.fig-dir "np4-schematic")))

(defn plot-basis-schematic [np &optional [annotate False]]
  (import islet)
  (sv c (futils.get-context)
      gll-best 1
      pats {0"--" gll-best "-"}
      clrs "krbcgm"
      x (npy.linspace -1 1 512)
      isl (islet.Islet)
      fs 12)
  (with [(pl-plot (, 5 (if annotate 3.33 2.5))
                  (+ c.fig-dir "basis-schematic-np" (str np)
                     (if annotate "-annotated" "")))]
    (for [method (, 0 1)]
      (sv y (.transpose (isl.eval method np x))
          pat (get pats method))
      (for [i (range np)]
        (pl.plot x (get y i) (+ (nth clrs (% i (len clrs))) pat)
                 :label (if (zero? i)
                            (get {0 "Natural GLL" gll-best "Islet GLL"}
                                 method))))
      (my-grid)
      (sv d 1.04)
      (pl.xlim (, (- d) d))
      (pl.ylim (, (if annotate -0.64 -0.22) 1.03)))
    (pl.xlabel "Reference-element coordinate" :fontsize fs)
    (pl.ylabel "Basis function value" :fontsize fs)
    (pl.text -1.22 1.1 (.format "Basis functions, $n_p$ = {}" np) :fontsize fs)
    (pl.figlegend :loc (, 0.44 0.9) :fontsize (- fs 1) :ncol 2)
    (when annotate
      (sv npa npy.array
          xgll (npa (poly.get-gll-x np))
          ireg 1
          xs (cut xgll ireg (+ ireg 2))
          xc (cc xs)
          clr "g" y -0.3 w 0.006 lw 2 ones (npa [1 1]))
      (for [i (range 2)]
        (pl.arrow (nth xgll (+ ireg i)) y 0 0.1 :width w :color clr))
      (pl.plot xs (* y ones) "-" (* xc ones) [-0.36 y] :color clr :lw lw)
      (pl.text xc -0.45 (.format "Region {}" ireg)
               :color clr :ha "center" :fontsize fs)
      (when (= np 6)
        (sv support [0 1 2 3 5]
            clr "r")
        (for [i support]
          (pl.arrow (nth xgll i) -0.5 0 0.1 :width w :color clr))
        (sv x (nth xgll 4))
        (pl.plot x -0.45 (+ clr "x") :markersize 14)
        (pl.text 0 -0.59 (.format "Support nodes for region {}" ireg)
                 :color clr :ha "center" :fontsize fs)))))

(when-inp ["basis-schematic" {:np int}]
  (plot-basis-schematic np :annotate (= np 6)))

(when-inp ["illustrations"]
  (sv c (futils.get-context))
  (illustrate-grids (+ c.fig-dir "illustrate-grids"))
  (draw-islet-timestep-schematic (+ c.fig-dir "islet-timestep-schematic")))

(when-inp ["islet-timestep-schematic"]
  (sv c (futils.get-context))
  (draw-islet-timestep-schematic (+ c.fig-dir "islet-timestep-schematic")))

(defn ms-draw-block [np x-beg y-beg ulen clr lw fs]
  (sv n (dec np)
      xl (* np ulen)
      yl (* n ulen)
      pat (+ clr "-"))
  (pl.plot [x-beg (+ x-beg (* np ulen))]
           [y-beg y-beg]
           pat :lw lw)
  (pl.plot [x-beg (+ x-beg (* np ulen))]
           [(+ y-beg yl) (+ y-beg yl)]
           pat :lw lw)
  (pl.plot [x-beg x-beg]
           [y-beg (+ y-beg yl)]
           pat :lw lw)
  (pl.plot [(+ x-beg xl) (+ x-beg xl)]
           [y-beg (+ y-beg yl)]
           pat :lw lw)
  (sv x (+ x-beg (- xl ulen)))
  (pl.plot [x x]
           [y-beg (+ y-beg yl)]
           (+ clr "--")  :lw lw)
  (pl.text (+ x-beg (* 0.5 n ulen)) (+ y-beg (* 0.5 yl))
           "$\\bar{B}$"
           :fontsize fs :ha "center" :va "center" :color clr)
  (pl.text (+ x-beg (* (+ n 0.5) ulen)) (+ y-beg (* 0.5 yl))
           "$b$"
           :fontsize fs :ha "center" :va "center" :color clr))

(defn ms-write-numbers [np down x-beg y-beg dx dy fs]
  (sv n (+ down (* 2 (dec np))))
  (for [i (range n)]
    (pl.text (+ x-beg (* i dx))
             (+ y-beg (* i dy))
             (str i)
             :fontsize fs :ha "center" :va "center")))

(defn ms-draw-elements [np x-beg y ulen line-pat node-pat fill-style]
  (sv n (dec np)
      e-len (* n ulen)
      x-len (* ulen (* 2 n))
      x-end (+ x-beg x-len)
      e-bdys [x-beg (+ x-beg (* n ulen)) x-end]
      xgll (poly.get-gll-x np))
  (pl.plot [x-beg x-end] [y y] line-pat
           [x-end (* (dec (* 2 np)) ulen)] [y y] (+ (first line-pat) ":")
           e-bdys [y y y]
           :lw 4 node-pat :fillstyle fill-style :markersize 16 :mew 4)
  (sv xs e-bdys)
  (for [i (range 2)]
    (sv xphys (cut (+ x-beg (* i e-len) (* (/ (+ 1 (npy.array xgll)) 2) e-len))
                   1 -1))
    (.extend xs (list xphys))
    (pl.plot xphys (* y (npy.ones (len xphys)))
             node-pat :fillstyle fill-style :markersize 11 :mew 4))
  (.sort xs)
  xs)

(defn ms-annotate-elements [xs y ylbl lbl fs]
  (for [i (range (len xs))]
    (pl.text (nth xs i) y (str i) :fontsize fs :ha "center" :va "center"))
  (pl.text (first xs) ylbl lbl :fontsize fs :ha "left" :va "center"))

(defn make-matrix-schematic [&optional np down filename]
  (sv c (futils.get-context))
  (svifn np 4 down 1 filename (+ c.fig-dir "matrix-schematic"))
  (assert (>= np 4))
  (assert (>= down 0))

  (sv n (dec np)
      ulen 1 hulen (/ ulen 2)
      L (* (+ down (* 2 (dec np))) ulen)
      fs 20 fsb (int (* fs 2))
      xgll (poly.get-gll-x np)
      xlate (- (first xgll) (nth xgll down)))
  
  (with [(pl-plot (, 6 7) filename)]
    (pl.plot [0 L] [L L] "k-"
             [0 0] [0 L] "k-")
    (ms-write-numbers np down hulen (+ L hulen) 1 0 fs)
    (ms-write-numbers np down (- hulen) (- L hulen) 0 -1 fs)
    (ms-draw-block np 0 (- L (* ulen (+ down n))) ulen "r" 2 fsb)
    (ms-draw-block np n 0 ulen "r" 2 fsb)
    (sv xs-src (ms-draw-elements np hulen (* -2 ulen) ulen
                                 "b-" "bo" "full")
        xs-tgt (ms-draw-elements np (+ hulen xlate) (* -1.35 ulen) ulen
                                 "g--" "go" "none"))
    (ms-annotate-elements xs-src (* -2.5 ulen) (* -3 ulen)
                          "Source (Eulerian)" fs)
    (ms-annotate-elements xs-tgt (* -0.85 ulen) (* -0.35 ulen)
                          "Target (Lagrangian)" fs)
    (pl.axis "off")))

(when-inp ["matrix-schematic"]
  (make-matrix-schematic))

(defn get-unit-cube-corners []
  (npy.array [[1 -1 -1] [1 1 -1] [-1 1 -1] [-1 -1 -1]
              [1 -1 1] [1 1 1] [-1 1 1] [-1 -1 1]]))

(defn normalize [v] (/ v (npy.linalg.norm v)))

(defn cs-get-cubedsphere-ne2-lines [n]
  (assert (>= n 2))
  (sv npa npy.array
      c (get-unit-cube-corners))
  (defn avg [&rest args]
    (sv p (.copy (nth c (first args))))
    (for [i (range 1 (len args))]
      (+= p (nth c (nth args i))))
    (/ p (len args)))
  (defn arc [a b]
    (npy.concatenate (, (col-vec (npy.linspace (nth a 0) (nth b 0) n))
                        (col-vec (npy.linspace (nth a 1) (nth b 1) n))
                        (col-vec (npy.linspace (nth a 2) (nth b 2) n)))
                     :axis 1))
  (defn add [sq d]
    (sv sqd (.copy sq))
    (for [e sqd] (+= e d))
    sqd)
  (defn xform [v i is j js k ks]
    (sv x [])
    (for [e v]
      (.append x (npy.concatenate
                  (, (* is (col-vec (getc e i)))
                     (* js (col-vec (getc e j)))
                     (* ks (col-vec (getc e k))))
                  :axis 1)))
    x)
  (sv sq (npy.concatenate (, (arc (nth c 0) (avg 0 1))
                             (cut (arc (avg 0 1) (avg 0 1 4 5)) 1)
                             (cut (arc (avg 0 1 4 5) (avg 0 4)) 1)
                             (cut (arc (avg 0 4) (nth c 0)) 1))
                          :axis 0))
  (sv side [] v [])
  (for [d (, (add sq [0 0 0])
             (add sq [0 1 0])
             (add sq [0 1 1])
             (add sq [0 0 1]))]
    (.append side d))
  (.extend v (xform side 0 1 1 1 2 1))
  (.extend v (xform side 1 -1 0 1 2 1))
  (.extend v (add (xform side 0 1 1 -1 2 1) [-2 0 0]))
  (.extend v (xform side 1 1 0 -1 2 1))
  (.extend v (xform side 1 1 2 1 0 1))
  (.extend v (add (xform side 1 1 2 1 0 1) [0 0 -2]))
  (for [e v] (normalize-rows! e))
  v)

(defn cs-make-separate-lines [lines]
  ;; make each arc a bunch of separate lines.
  (sv sl [])
  (for [e lines]
    (for [i (range (dec (len e)))]
      (.append sl (get e (, [i (inc i)] (s-all))))))
  sl)

(defn cs-project-lines [nml lines]
  (sv vfront [] vback []
      xhat (normalize (npy.cross [0 0 1] nml))
      yhat (col-vec (normalize (npy.cross nml xhat)))
      xhat (col-vec xhat)
      nml (col-vec (normalize nml)))
  (for [e lines]
    (sv p (npy.concatenate (, (npy.dot e xhat)
                              (npy.dot e yhat))
                           :axis 1)
        s (npy.dot e nml)
        fmask (vectorize (>= s 0))
        bmask (npy.logical-not fmask))
    (.append vfront (get p (, fmask (s-all))))
    (.append vback  (get p (, bmask (s-all)))))
  (, vfront vback xhat yhat nml))

(defn cs-ref-to-sphere [crnrs a b]
  (sv oma (- 1 a)
      omb (- 1 b))
  (normalize (+ (* a b (nth crnrs 0))
                (* oma b (nth crnrs 1))
                (* oma omb (nth crnrs 2))
                (* a omb (nth crnrs 3)))))

(defn make-gll-pts [crnrs np]
  (sv r (poly.get-gll-x np)
      ps [])
  (for [i (range np)]
    (sv a (/ (+ 1 (nth r i)) 2))
    (for [j (range np)]
      (sv b (/ (+ 1 (nth r j)) 2))
      (.append ps (cs-ref-to-sphere crnrs a b))))
  ps)

(defn cs-make-physgrid [crnrs nf n]
  (sv lns []
      zon (npy.linspace 0 1 n))
 (for [i (range (inc nf))]
    (sv a (/ i nf)
        ln (npy.zeros (, n 3)))
    (for [(, j b) (enumerate zon)]
      (sv (get ln j) (cs-ref-to-sphere crnrs a b)))
    (.append lns ln))
 (for [i (range (inc nf))]
    (sv b (/ i nf)
        ln (npy.zeros (, n 3)))
    (for [(, j a) (enumerate zon)]
      (sv (get ln j) (cs-ref-to-sphere crnrs a b)))
    (.append lns ln))
  lns)

(defn cubed-sphere-subelem-grid-schematic [&optional npv npt filename nf]
  (sv c (futils.get-context))
  (svifn npv 4 npt 6
         filename (+ c.fig-dir "cubed-sphere-subelem-grid-schematic"))
  (assert (and (>= npv 4) (>= npt npv)))

  (sv vlines (cs-get-cubedsphere-ne2-lines 50)
      vlines (cs-make-separate-lines vlines)
      nml-ref [0.3 -0.7 0.2]
      (, vfront vback xhat yhat nml) (cs-project-lines nml-ref vlines))

  (defn project-pts [g]
    (sv p (npy.zeros (, (len g) 2)))
    (for [(, i e) (enumerate g)]
      (sv (get p (, i (s-all))) [(npy.float64 (npy.dot e (col-vec xhat)))
                                 (npy.float64 (npy.dot e (col-vec yhat)))]))
    p)
  (sv crnrs (lfor e [[1 -1 1] [0 -1 1] [0 -1 0] [1 -1 0]]
                  (row-vec (normalize e)))
      gv (project-pts (make-gll-pts crnrs npv))
      gt (project-pts (make-gll-pts crnrs npt)))

  (with [(pl-plot (, 6 6) filename)]
    (sv ax (pl.gca)
        lim (, -1.1 1.1)
        g 0.8
        lcb (mpl.collections.LineCollection vback  :color [g g g] :lw 2)
        lcf (mpl.collections.LineCollection vfront :color [0 0 0] :lw 2)
        theta (npy.linspace 0 (* 2 math.pi) 200))
    (ax.add-collection lcb)
    (sv g 0.4)
    (pl.plot (npy.cos theta) (npy.sin theta) "--" :lw 2
             :color (if (none? nf) "g" [g g g]))
    (ax.add-collection lcf)
    (unless (none? nf)
      (sv pgls (cs-make-separate-lines (cs-make-physgrid crnrs nf 50))
          (, pgls - - - -) (cs-project-lines nml-ref pgls)
          lpg (mpl.collections.LineCollection pgls :color "g" :lw 1))
      (ax.add-collection lpg))
    (pl.plot (getc gv 0) (getc gv 1) "bo" :markersize 12)
    (pl.plot (getc gt 0) (getc gt 1) "ro" :markersize 6)
    (pl.xlim lim) (pl.ylim lim)
    (pl.axis "square")
    (pl.axis "off")))

(when-inp ["cubed-sphere-subelem-grid-schematic"]
  (cubed-sphere-subelem-grid-schematic :nf 2))

(defn make-gll-grid [xb np &optional [end False]]
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

(defn ones [x] (npy.ones x.shape))

(defn wrap-01 [xi]
  (sv mask (< xi 0)
      (get xi mask) (+ (get xi mask) 1)
      mask (> xi 1)
      (get xi mask) (- (get xi mask) 1))
  xi)

(defn interp-sup [x y isup xi]
  (sv xsup (npy.array (lfor i isup (get x i)))
      ysup (npy.array (lfor i isup (get y i))))
  (poly.eval-lagrange-poly xsup ysup xi))

(defn csl-cubic-interp [x y xi]
  (sv nc (dec (len x))
      yi (npy.zeros xi.shape)
      xi (wrap-01 xi)
      nodes [])
  (for [ic (range nc)]
    (sv mask (npy.logical-and (>= xi (get x ic)) (< xi (get x (inc ic)))))
    (.append nodes (list (first (npy.nonzero mask))))
    (unless (npy.any mask) (continue))
    (sv isup [(% (dec ic) nc) ic (inc ic) (% (+ ic 2) nc)]
        (get yi mask) (interp-sup x y isup (get xi mask))))
  (, yi nodes))

(defn gll-interp [np x y xi]
  (sv d (dec np)
      nc (// (dec (len x)) d)
      yi (npy.zeros xi.shape)
      xi (wrap-01 xi)
      nodes [])
  (for [ic (range nc)]
    (sv xlo (get x (* ic d))
        xhi (get x (* (inc ic) d))
        mask (npy.logical-and (>= xi xlo) (< xi xhi)))
    (.append nodes (list (first (npy.nonzero mask))))
    (unless (npy.any mask) (continue))
    (sv isup (list (range (* ic d) (inc (* (inc ic) d))))
        (get yi mask) (interp-sup x y isup (get xi mask))))
  (, yi nodes))

(defn isl-1d-schematic [&optional ne-gll ne-csl np filename]
  (sv c (futils.get-context))
  (svifn ne-gll 5 np 4
         filename (+ c.fig-dir "isl-1d-schematic"))
  (svifn ne-csl (* ne-gll (dec np)))
  (assert (>= np 4))

  (defn f-fn [x]
    (sv ctr 0.4
        fac 3
        f (* 0.5 (+ 1 (npy.cos (* math.pi (* fac (- x ctr))))))
        mask (npy.logical-or (< x (- ctr (/ fac))) (> x (+ ctr (/ fac))))
        (get f mask) 0)
    f)

  (sv csl-grid (npy.linspace 0 1 (inc ne-csl))
      el-grid (npy.linspace 0 1 (inc ne-gll))
      gll-grid (make-gll-grid el-grid np :end True))

  (defn draw-csl-grid [xos yos]
    (pl.plot (+ xos (npy.array [0 1])) [yos yos] "k-")
    (pl.plot (+ xos csl-grid) (* yos (ones csl-grid)) "ko"
             :markersize 4 :fillstyle "full"))

  (defn draw-gll-grid [xos yos]
    (pl.plot (+ xos el-grid) (* yos (ones el-grid)) "ko"
             :markersize 7 :fillstyle "full")
    (pl.plot (+ xos (npy.array [0 1])) [yos yos] "k-")
    (pl.plot (+ xos gll-grid) (* yos (ones gll-grid)) "ko"
             :markersize 4 :fillstyle "full"))

  (sv xlate -0.235
      f-csl (f-fn csl-grid)
      (, i-csl nodes-csl) (csl-cubic-interp csl-grid f-csl (+ csl-grid xlate))
      f-gll (f-fn gll-grid)
      (, i-gll nodes-gll) (gll-interp np gll-grid f-gll (+ gll-grid xlate)))
  (print nodes-csl)

  (with [(pl-plot (, 8 3.5) filename)]
    (sv gap 0.1
        yos 1.5
        tfs 12
        p 0.5
        gray [p p p])

    (pl.text 0.5 (+ yos 1.2) "Classical ISL"
             :fontsize tfs :ha "center")
    (pl.text (+ 1.5 gap) (+ yos 1.2) "Element-based ISL"
             :fontsize tfs :ha "center")

    (pl.text 0 (+ yos 0.9) "Time step n+1" :fontsize tfs)
    (pl.text 0 0.9 "Time step n" :fontsize tfs)
    (for [i (get nodes-csl 5)]
      (pl.plot [(+ (get csl-grid i) xlate) (get csl-grid i)]
               [0 yos] "g--")
      (pl.plot (+ (get csl-grid i) xlate) 0 "gx")
      (pl.plot (get csl-grid i) yos "gs" :fillstyle "none"))

    (sv i 9
        d -0.2)
    (pl.plot (get csl-grid i) (+ yos (get i-csl i)) "ro"
             :fillstyle "none")
    (pl.plot [(+ (get csl-grid i) xlate) (get csl-grid i)]
             [d d] "g-")
    (pl.plot (get csl-grid i) d "gs" :fillstyle "none")
    (pl.plot (+ (get csl-grid i) xlate) d "gx"
             :fillstyle "none")
    (pl.text (+ (get csl-grid i) 0.02) d
             "Translation distance $\Delta x$"
             :fontsize tfs :va "center")

    (draw-csl-grid 0 0)
    (sv xi (npy.linspace 0 1 111)
        (, yi -) (csl-cubic-interp csl-grid f-csl xi))
    (pl.plot xi yi "-" :color gray)
    (sv isup (list (range 4 8))
        xi (npy.linspace (get csl-grid (first isup))
                         (get csl-grid (last isup))
                         30)
        yi (interp-sup csl-grid f-csl isup xi))
    (pl.plot xi yi "r:" :lw 2)
    (sv xi (npy.linspace (get csl-grid (get isup 1))
                         (get csl-grid (get isup 2))
                         10)
        yi (interp-sup csl-grid f-csl isup xi))
    (pl.plot xi yi "r-" :lw 2)
    (pl.plot (get csl-grid isup) (get f-csl isup) "ro"
             :fillstyle "none")
    (pl.plot csl-grid f-csl "k.")
    
    (draw-csl-grid 0 yos)
    (pl.plot csl-grid (+ yos i-csl) "k.")

    (sv ni 1 dx (+ 1 gap))
    (draw-gll-grid (+ dx) 0)
    (draw-gll-grid (+ dx) yos)

    (for [i (get nodes-gll ni)]
      (pl.plot [(+ (get gll-grid i) xlate dx) (+ (get gll-grid i) dx)]
               [0 yos] "g--")
      (pl.plot (+ (get gll-grid i) xlate dx) 0 "gx")
      (pl.plot (+ (get gll-grid i) dx) yos "gs" :fillstyle "none"))
    
    (sv xi (npy.linspace 0 1 111)
        (, yi -) (gll-interp np gll-grid f-gll xi))
    (pl.plot (+ dx xi) yi "-" :color gray)
    (sv isup (list (range (* (dec np) ni) (inc (* (dec np) (inc ni)))))
        xi (npy.linspace (get gll-grid (first isup))
                         (get gll-grid (last isup))
                         30)
        yi (interp-sup gll-grid f-gll isup xi))
    (pl.plot (+ dx xi) yi "b-" :lw 2)
    (pl.plot (+ dx (get gll-grid isup)) (get f-gll isup) "bo"
             :fillstyle "none")
    (pl.plot (+ dx gll-grid) f-gll "k.")
    
    (pl.plot (+ dx gll-grid) (+ yos i-gll) "k.")
    (sv idxs (get nodes-gll ni))
    (pl.plot (+ dx (get gll-grid idxs)) (+ yos (get i-gll idxs)) "bo"
             :fillstyle "none")
    
    (pl.axis "off")))

(when-inp ["isl-1d-schematic"]
  (isl-1d-schematic))
