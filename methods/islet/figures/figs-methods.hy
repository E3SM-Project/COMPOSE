(require [amb3 [*]])
(import amb3 [amb3 [*]] [figsutils :as futils]
        [scipy.linalg :as linalg]
        math re poly)

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
    (write "\\begin{tabular}{r|c|l|l}\n")
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

(defn illustrate-grids [img-fname]
  (defn draw-gll-dots [np marker]
    (sv xgll (col-vec (poly.get-gll-x np))
        e (npy.ones (, np 1))
        X (npy.dot e (.transpose xgll))
        Y (npy.dot xgll (.transpose e)))
    (pl.plot X Y marker :markersize 12))
  (sv nf 6 elc "k")
  (with [(pl-plot (, 5 5) img-fname :format "pdf")]
    (pl.plot [-1 1] [1 1] elc [-1 1] [-1 -1] elc [1 1] [-1 1] elc [-1 -1] [-1 1] elc
               :linewidth 2
               :color (if 0 "#1D4BFA" "#3838FF"))
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
    (draw-gll-dots 8 "r.")
    (pl.axis "equal")
    (pl.axis "off")))

(defn draw-np4-schematic [img-fname]
  (import islet)
  (sv np 4 fs 12
      x-gll (poly.get-gll-x np)
      nx 256
      x (npy.linspace -1 1 nx)
      isl (islet.Islet)
      (, yn yo yb) (lfor m (, 0 1 3) (.transpose (isl.eval m np x)))
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
    (pl.xlabel "Reference coordinate" :fontsize fs)))

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

(defn plot-slmmir-vs-heuristic [c d img-fname
                                &optional nps prop-preserve ic norm pum-thr lebesgue]
  (svifn nps (list (range 4 14)) prop-preserve False ic "gau" norm :l2
         pum-thr 1e-6 lebesgue False)
  (sv npa npy.array fs 11
      plot (if lebesgue pl.semilogy pl.loglog))
  (print img-fname)
  (with [(pl-plot (, 4 4.2) img-fname :format "pdf")]
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
            (pl.xlabel (+ "$\log_{10}$ " (norm->str norm) " heuristic") :fontsize fs)))
    (my-grid)
    (pl.title (.format (+ "Test problem vs.~heuristic:\n"
                          "nondivergent flow, 1.5$^\circ$, {}, long steps\n"
                          "$p$-refinement, {}")
                       (futils.ic-short2long ic)
                       (+ (if prop-preserve "" "no ") "property preservation"))
              :fontsize fs)
    (pl.legend :loc "best" :fontsize fs)
    (pl.axis "tight")
    (pl.ylabel (+ "$\log_{10}$ " (norm->str norm) " relative error") :fontsize fs)))

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
    (make-reference-slope-triangle [(/ 2.8e-3 f) (* 2.8e-2 f)]
                                   [(* 3e-8 (** f slope)) (/ 3e-12 (** f slope))]
                                   slope "k-"
                                   :opposite True
                                   :kwargs-text {"fontsize" (+ fs 4)})
    (set-loglog-ticks)
    (my-grid)
    (pl.legend :loc "upper left" :fontsize (dec fs) :ncol 1)
    (pl.xlabel "$\log_{10}$ Element size relative random perturbation $\delta$"
               :fontsize fs)
    (pl.ylabel "$\log_{10}$ (max $|\lambda|$ - 1)" :fontsize fs)
    (pl.title "Perturbed uniform mesh metric" :fontsize fs)
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
            (nth ms 2) "GLL nodal subset"})
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
    (pl.xlabel "Translation, $\Delta x$, relative to element size 1" :fontsize fs)
    (pl.title (.format "$n_p$ {} methods" np) :fontsize fs)
    (pl.legend :loc (, 0.02 0.15) :fontsize fs)
    (pl.xticks (npy.linspace 0 1 11))
    (my-grid)))

;;; drivers

(when-inp ["dev-parse"]
  (sv fname "data/search-0.txt"
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
    (.extend blns (parse-search-list (+ "data/" fname))))
  (sv blns (uniquify-search-list blns))
  (write-slmmir-script blns script-fname))

(when-inp ["plot-slmmir-vs-heuristic"]
  (sv fnames (, "slmmir-on-basis-lines-2.txt")
      c (futils.get-context)
      lebesgue False
      d {})
  (for [fname fnames]
    (sv d (parse-slmmir-output (+ "data/" fname) :d d :lebesgue lebesgue)))
  (for [(, norm pp ic) (, (, :l2 False "gau") (, :l2 True "cos"))]
    (plot-slmmir-vs-heuristic
      c d (.format (+ "{}slmmir-vs-heuristic-{}-{}-{}" (if lebesgue "-leb" ""))
                   c.fig-dir ic (if pp "pp" "nopp") (name norm))
      :nps [6 7 8 9 10]
      :prop-preserve pp :norm norm :ic ic :lebesgue lebesgue)))

(when-inp ["illustrations"]
  (sv c (futils.get-context))
  (illustrate-grids (+ c.fig-dir "illustrate-grids")))

(when-inp ["tables"]
  (for [xnodes (, 'gll)]
    (write-methods-latex-tables (parse-cpp-methods "islet-methods.txt")
                                xnodes (.format "figs/methods-table-{}.tex"
                                                (name xnodes)))))

(when-inp ["pum-vs-perturb"]
  (sv fname "pum_perturb_plot-041021.txt"
      d (parse-pum-vs-perturb (+ "data/" fname))
      c (futils.get-context))
  (plot-pum-vs-perturb c d (+ c.fig-dir "pum-vs-perturb")))

(when-inp ["meam1-and-pum-vs-dx"]
  (sv c (futils.get-context)
      data-dir "data/"
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
      pats {0 "--" 3 "-"}
      clrs "krbcgm"
      x (npy.linspace -1 1 512)
      isl (islet.Islet)
      fs 12)
  (with [(pl-plot (, 6 (if annotate 4 3))
                  (+ c.fig-dir "basis-schematic-np" (str np)
                     (if annotate "-annotated" "")))]
    (for [method (, 0 3)]
      (sv y (.transpose (isl.eval method np x))
          pat (get pats method))
      (for [i (range np)]
        (pl.plot x (get y i) (+ (nth clrs (% i (len clrs))) pat)
                 :label (if (zero? i) (get {0 "Natural GLL" 3 "Islet GLL"} method))))
      (my-grid)
      (sv d 1.04)
      (pl.xlim (, (- d) d))
      (pl.ylim (, (if annotate -0.64 -0.22) 1.03)))
    (pl.xlabel "Reference coordinate" :fontsize fs)
    (pl.ylabel "Basis function value" :fontsize fs)
    (pl.text -1 1.1 (.format "Basis functions, $n_p$ = {}" np) :fontsize fs)
    (pl.figlegend :loc (, 0.49 (if annotate 0.91 0.88)) :fontsize fs :ncol 2)
    (when annotate
      (sv npa npy.array
          xgll (npa (poly.get-gll-x np))
          ireg 1
          xs (cut xgll ireg (+ ireg 2))
          xc (cc xs)
          clr "g" y -0.3 w 0.006 lw 2 ones (npa [1 1]))
      (for [i (range 2)]
        (pl.arrow (nth xgll (+ ireg i)) y 0 0.1 :width w :color clr))
      (pl.plot xs (* y ones) (+ clr "-") (* xc ones) [-0.36 y] :color clr :lw lw)
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
  (plot-basis-schematic np :annotate True))
