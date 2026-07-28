// MathJax configuration for `pymdownx.arithmatex` in generic mode.
//
// Arithmatex wraps math in `<script type="math/tex">`-free `\(...\)` and
// `\[...\]` delimiters, so those are what MathJax is told to look for. The
// `subscriptionend`/`document$` hook re-typesets after Material's instant
// navigation swaps page content without a reload.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
