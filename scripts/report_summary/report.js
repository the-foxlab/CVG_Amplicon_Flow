// One global sample selector drives every section at once: switching it shows
// only the panels of the selected sample across the whole report. With a single
// sample no selector is rendered and all panels stay visible.
const select = document.querySelector('[data-sample-select]');
if (select) {
	const update = () => {
		document.querySelectorAll('.sample-panel').forEach((panel) => {
			panel.hidden = panel.dataset.sample !== select.value;
		});
	};
	select.addEventListener('change', update);
	update();
}