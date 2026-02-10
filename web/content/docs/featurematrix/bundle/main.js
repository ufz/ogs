import {
  buildFeatureTree,
  buildFileTree,
  filterData,
  collectHighlightedLines,
} from "./utils.js";

let data = [];
let tree = {};
let fileTree = {};
const fileColors = {};
let selectedFiles = [];
const palette = [
  "red",
  "blue",
  "green",
  "orange",
  "purple",
  "teal",
  "magenta",
  "brown",
  "navy",
  "lime",
];

async function init() {
  console.log("init: calling");
  try {
    const resp = await fetch("features.json");
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    return await resp.json();
  } catch (err) {
    console.error("Failed to load features.json:", err);
    alert("Could not load features.json. Check console for details.");
    return [];
  }
}

document.addEventListener("DOMContentLoaded", async () => {
  data = await init();
  tree = buildFeatureTree(data);
  fileTree = buildFileTree(data);
  renderFeatures();
  update();
});

function createFeatureCheckbox(group, subfeature, onChange, onHover, onLeave) {
  const label = document.createElement("label");
  const cb = document.createElement("input");
  cb.type = "checkbox";
  cb.value = `${group}:${subfeature}`;
  cb.addEventListener("change", onChange);
  const name = document.createTextNode(subfeature);
  const dots = document.createElement("div");
  dots.className = "dots-container";
  const count = document.createElement("span");
  count.className = "count";
  count.textContent = "0";
  label.append(cb, name, dots, count);
  label.addEventListener("mouseenter", () => onHover(cb.value));
  label.addEventListener("mouseleave", onLeave);
  return label;
}

function createFeatureGroup(groupName, subfeatures, checkboxFactory) {
  const det = document.createElement("details");
  det.open = true;
  const sum = document.createElement("summary");
  sum.textContent = groupName;
  det.append(sum);
  Array.from(subfeatures)
    .sort()
    .forEach((s) => det.append(checkboxFactory(groupName, s)));
  return det;
}

function renderFeatures() {
  console.log("rendering features:", tree);
  const container = document.getElementById("features");
  container.innerHTML = "";
  Object.entries(tree)
    .sort()
    .forEach(([group, subfeatures]) => {
      const groupElem = createFeatureGroup(group, subfeatures, (g, s) =>
        createFeatureCheckbox(
          g,
          s,
          () => update(),
          (f) => hoverFeature(f),
          () => clearHover(),
        ),
      );
      container.append(groupElem);
    });
}

function getSelectedFilters() {
  return Array.from(document.querySelectorAll("#features input:checked")).map(
    (i) => i.value,
  );
}

function getAvailableFeatures(data) {
  return new Set(data.flatMap((d) => d.features));
}

function update() {
  console.log("updating from data:", data.length);
  const filters = getSelectedFilters();
  const filtered = filterData(data, filters);

  // Prune selectedFiles to visible
  selectedFiles = selectedFiles.filter((f) =>
    filtered.some((d) => d.file === f),
  );
  renderFiles(filtered);
  // Show only available features and update counts
  const available = getAvailableFeatures(filtered);
  document.querySelectorAll("#features details").forEach((det) => {
    let anyVisible = false;
    det.querySelectorAll("label").forEach((label) => {
      const cb = label.querySelector("input");
      const f = cb.value;
      const cnt = filtered.filter((d) => d.features.includes(f)).length;
      label.querySelector(".count").textContent = cnt;
      if (available.has(f)) {
        label.style.display = "flex";
      } else {
        if (cb.checked) cb.checked = false;
        label.style.display = "none";
      }
      if (label.style.display !== "none") anyVisible = true;
    });
    det.style.display = anyVisible ? "" : "none";
  });
  clearDots();
  emphasize();
}

function getFileColor(file) {
  if (!fileColors[file]) {
    const index = Object.keys(fileColors).length;
    fileColors[file] = palette[index % palette.length];
  }
  return fileColors[file];
}

function createFileEntry(
  fn,
  featureCount,
  featureCoverage,
  color,
  selected,
  onClick,
  onHover,
  onLeave,
) {
  if (featureCoverage === undefined) {
    throw new Error(`Missing feature_coverage for file: ${fn}`);
  }
  const div = document.createElement("div");
  div.className = "file-item";
  div.dataset.file = fn;

  if (selected) div.classList.add("selected");

  const countSpan = document.createElement("span");
  countSpan.className = "file-count";
  countSpan.textContent = featureCount;

  const fcSpan = document.createElement("span");
  fcSpan.className = "file-coverage";
  fcSpan.textContent = featureCoverage;

  const basename = fn.includes("/")
    ? fn.substring(fn.lastIndexOf("/") + 1)
    : fn;
  const nameNode = document.createTextNode(basename);

  // Order: feature count | feature coverage | file name
  div.append(countSpan, fcSpan, nameNode);

  div.style.setProperty("--file-color", color);
  div.addEventListener("click", () => onClick(div, fn));
  div.addEventListener("mouseenter", () => onHover(fn));
  div.addEventListener("mouseleave", onLeave);

  return div;
}

function renderFiles(list) {
  const filesDiv = document.getElementById("files");
  const availableFeatures = getAvailableFeatures(list);
  filesDiv.innerHTML = "";

  Object.keys(fileTree)
    .sort()
    .forEach((path) => {
      const groupName = path === "." ? "." : path;
      const files = fileTree[path].filter((f) =>
        list.some((d) => d.file === f),
      );
      if (!files.length) return;

      const itemFactory = (_group, file) => {
        const entry = data.find((d) => d.file === file);
        const featureCount = entry.features.filter((f) =>
          availableFeatures.has(f),
        ).length;
        const color = getFileColor(file);
        const selected = selectedFiles.includes(file);
        const featureCoverage = Number(entry.feature_coverage).toFixed(2);
        return createFileEntry(
          file,
          featureCount,
          featureCoverage,
          color,
          selected,
          toggleSelect,
          hoverFile,
          clearHover,
        );
      };

      const group = createFeatureGroup(groupName, files, itemFactory);
      filesDiv.append(group);
    });

  if (!filesDiv.children.length) filesDiv.innerHTML = "<div>No matches</div>";
}

async function showProjectFile(file, lines) {
  const viewer = document.getElementById("project-viewer");
  viewer.innerHTML = `<div class="project-filename">${file}</div>`;
  try {
    const resp = await fetch(`Tests/Data/${file}`);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    const text = await resp.text();

    // Prism highlight
    const highlighted = Prism.highlight(text, Prism.languages.xml, "xml");
    const highlightedLines = highlighted.split("\n");

    const greenLines = collectHighlightedLines(lines);

    let html = "";
    for (let i = 0; i < highlightedLines.length; ++i) {
      const ln = i + 1;
      const bar = greenLines.has(ln) ? `<span class="green-bar"></span>` : "";
      // Use <span style="white-space: pre"> to preserve spaces in each line
      html += `<div class="file-line">${bar}<span style="white-space:pre">${highlightedLines[i] || "\u00A0"}</span></div>`;
    }

    viewer.innerHTML = `
      <div class="project-filename">${file}</div>
      <pre class="file-pre"><code class="language-xml">${html}</code></pre>
    `;
  } catch (err) {
    viewer.innerHTML = `<div class="project-filename">${file}</div><pre>Could not load file: ${file}\n${err}</pre>`;
  }
}

function toggleSelect(div, file) {
  const idx = selectedFiles.indexOf(file);
  if (idx > -1) {
    selectedFiles.splice(idx, 1);
    div.classList.remove("selected");
  } else {
    selectedFiles.push(file);
    div.classList.add("selected");
  }
  emphasize();
  if (selectedFiles.length) {
    const lastFile = selectedFiles[selectedFiles.length - 1];
    const entry = data.find((d) => d.file === lastFile);
    showProjectFile(file, entry.lines);
  } else {
    document.getElementById("project-viewer").textContent = "";
  }
}

function clearDots() {
  document
    .querySelectorAll(".dots-container")
    .forEach((dc) => (dc.innerHTML = ""));
}

function emphasize() {
  clearDots();
  document
    .querySelectorAll(".dots-container")
    .forEach(
      (dc) =>
        (dc.style.gridTemplateColumns = `repeat(${selectedFiles.length}, auto)`),
    );
  selectedFiles.forEach((file, idx) => {
    data
      .find((d) => d.file === file)
      .features.forEach((f) => {
        const inp = document.querySelector(`input[value="${CSS.escape(f)}"]`);
        const dc = inp.closest("label").querySelector(".dots-container");
        const dot = document.createElement("span");
        dot.className = "dot";
        dot.style.background = fileColors[file];
        dot.style.gridColumnStart = idx + 1;
        dc.append(dot);
      });
  });
}

function hoverFeature(f) {
  document.querySelectorAll(".file-item").forEach((div) => {
    const file = div.dataset.file;
    const d = data.find((d) => d.file === file);
    if (d?.features.includes(f)) div.classList.add("hovered");
  });
}

function hoverFile(file) {
  data
    .find((d) => d.file === file)
    .features.forEach((f) => {
      document
        .querySelector(`input[value="${CSS.escape(f)}"]`)
        .closest("label")
        .classList.add("hovered");
    });
}

function clearHover() {
  document
    .querySelectorAll(".file-item, label")
    .forEach((el) => el.classList.remove("hovered"));
}
