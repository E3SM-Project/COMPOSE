(require [amb3 [*]])
(import [amb3 [*]]
        [figsutils [*]]
        math colorsys)

(assoc matplotlib.rcParams "savefig.dpi" 300)
(do (pl-require-type1-fonts))

(defn get-errs [d ne np cycle]
  (get d "exact" "movingvortices" 1 "pcsl" "none" "none" 0 np "vor" cycle ne :C))

(defn make-vortex-cmap []
  (sv lsc (fn [colors lo hi n]
            ((matplotlib.colors.LinearSegmentedColormap.from-list
              "filament" colors 11)
             (npy.linspace lo hi n)))
      h2r colorsys.hsv-to-rgb
      s 0.2
      red (lsc [(h2r 0 1 1) (h2r 0 s 1)] 0 1 11)
      b (/ 2 3)
      blue (lsc [(h2r b 1 1) (h2r b s 1)] 1 0 11)
      cmap (matplotlib.colors.ListedColormap
            (npy.concatenate (, blue red) :axis 0)))
  cmap)

(defn fig-vortex [c d configs]
  (print "TODO config labels")
  (defn draw-circle []
    (sv th (npy.linspace 0 (* 2 math.pi) 512)
        r 0.995
        x (* r (npy.cos th)) y (* r (npy.sin th))
        g 0.8)
    (pl.plot x y "-" :lw 0.75 :color [g g g]))
  (sv ncol (len configs)
      fs 11
      cmap (make-vortex-cmap)
      ecmap (matplotlib.colors.LinearSegmentedColormap.from-list
             "filament" [[1 1 1] [0 0 0]]))
  (with [(pl-plot (, (* 2 ncol) 4) (+ c.fig-dir "vortex")
                  :format "png"
                  :tight False)]
    (for [(, icon config) (enumerate configs)]
      (sv (, ne np nstep cycle recon) config
          e (get-errs d ne np cycle)
          arrs (read-slmmir-io-arrays
                (+ c.data-dir f"vortex-imgs/ne{ne}np{np}nstep{nstep}{recon}.bin")))
      (sv img (get arrs (+ 1 (* 3 cycle)))
          diff (/ (npy.abs (- img (get arrs (+ 2 (* 3 cycle)))))
                  (npy.max (npy.abs (get arrs (+ 2 (* 3 cycle))))))
          imgs [img diff])
      (for [i (range 2)]
        (sv img (get imgs i)
            mask (!= img 0)
            ext (, (npy.min (get img mask)) (npy.max img)))
        (when (zero? i)
          (sv nmask (npy.logical-not mask)))
        (sv (get img nmask) npy.nan)
        (print "corner" (get img (, 0 0)))
        (print "ext" ext)
        (if (zero? i)
          (sv vmin 0.45 vmax 1.55)
          (sv vmax (get e :li)
              vmin 0))
        (sv rect [(/ icon ncol) (- 1 (* i 0.5))
                  (/ 0.97 ncol) (* 0.97 0.5)]
            ax (pl.axes rect))
        (pl.imshow img :cmap (if (zero? i) cmap ecmap)
                   :vmin vmin :vmax vmax
                   :origin "upper"
                   :extent [-1 1 -1 1])
        (pl.axis "off")
        (pl.text -1 -1
                 (.format "({})" (if (zero? i)
                                   (get "abcd" icon)
                                   (get "efgh" icon)))
                 :ha "left" :va "bottom" :fontsize fs)
        (draw-circle)
        (unless (zero? i)
          (pl.text 0 0.5
                   (.format "$n_e \, {}$, $n_p \, {}$, cycle {}, $n_{{step}} \, {}$\n"
                            ne np cycle (* cycle nstep))
                   :ha "center" :fontsize fs)
          (unless (= recon "constant")
            (for [k (range 3)]
              (pl.text -0.4 (+ 0.0 (* 0.2 (- 2 k)))
                       (.format (+ (get ["$l_1$" "$l_2$" "$l_{{\infty}}$"] k)
                                   " {:1.2e}")
                                (get e (get (, :l1 :l2 :li) k)))
                       :ha "right" :fontsize fs)))
          (when (= recon "constant")
            (pl.text 0 0.9 "constant in subcell" :ha "center" :fontsize fs)))
        (cond [(and (zero? i) (= icon 0))
               (sv f (pl.gcf)
                   cax (f.add-axes [(- (get rect 0) 0.04)
                                    (+ (get rect 1) (* 0.1 (get rect 3)))
                                    (* 0.1 (get rect 2))
                                    (* 0.8 (get rect 3))]))
               (pl.colorbar :cax cax :ticks [0.5 0.75 1 1.25 1.5]
                            :aspect 10)
               (cax.yaxis.set-ticks-position "left")]
              [(= i 1)
               (sv f (pl.gcf)
                   cax (f.add-axes [(+ (get rect 0) (* 0.1 (get rect 2)))
                                    (+ (get rect 1) 0.1)
                                    (* 0.8 (get rect 2))
                                    (* 0.05 (get rect 3))])
                   cb (pl.colorbar :cax cax :orientation "horizontal"
                                   :ticks [0 vmax]
                                   :shrink 0.6 :aspect 10))
               (cb.ax.set-xticklabels ["0" (.format "{:1.2e}" vmax)])])))))

(when-inp ["fig-vortex"]
  (sv c (jcp-context (get-context))
      d (acc-parse (+ c.data-dir "vortex-imgs.txt")
                   :map-nstepfac jcp-nenp2nstepfac)
      configs [(, 3 13 72 1 "constant")
               (, 3 13 72 1 "bilin")
               (, 6 13 144 1 "bilin")
               (, 12 13 288 2 "bilin")])
  (fig-vortex c d configs))
