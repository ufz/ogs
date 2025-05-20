$(document).ready(function () {
  MicroModal.init();

  $('#nav-link-search').click(function (ev) {
    ev.preventDefault();

    MicroModal.show('modal-2', {
      onClose: function () { $('.nav-link-contact').blur(); },
      disableFocus: true
    });

    document.querySelector('.pagefind-ui__search-input').focus();
  });

  // Trigger search popup for f without modifier keys only
  document.addEventListener('keydown', (event) => {
    if (event.key == 'f' &&
      !event.ctrlKey &&
      !event.altKey &&
      !event.metaKey &&
      !event.shiftKey) {

      const modals = document.querySelectorAll('#modal-2');
      modals.forEach(element => {
        if (!element.classList.contains('is-open')) {
          event.preventDefault();
        }
      });

      MicroModal.show('modal-2', {
        onClose: function () { $('.nav-link-contact').blur(); },
        disableFocus: false
      });
    }
  }, false);
});
