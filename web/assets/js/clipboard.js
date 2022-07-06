var clipboard = new ClipboardJS('.copy-code-btn');

function resetBtn(element) {
  element.innerHTML = "Copy";
 }

clipboard.on('success', function(e) {
  e.trigger.innerHTML = "Copied!";
  window.setTimeout(resetBtn, 1200, e.trigger)
  e.clearSelection();
});
