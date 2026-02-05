export function buildFeatureTree(data) {
  const tree = {};
  const featuresList = Array.from(new Set(data.flatMap((d) => d.features)));
  featuresList.forEach((f) => {
    const [g, s] = f.split(":");
    (tree[g] = tree[g] || new Set()).add(s);
  });
  // Convert Sets to Arrays for testing
  Object.keys(tree).forEach((g) => (tree[g] = Array.from(tree[g])));
  return tree;
}

export function buildFileTree(data) {
  const fileTree = {};
  data.forEach((d) => {
    const path = d.file.includes("/")
      ? d.file.substring(0, d.file.lastIndexOf("/"))
      : ".";
    (fileTree[path] = fileTree[path] || []).push(d.file);
  });
  return fileTree;
}

export function filterData(data, filters) {
  return data.filter((d) => filters.every((f) => d.features.includes(f)));
}

export function collectHighlightedLines(lines) {
  const set = new Set();
  for (const featureIntervals of lines) {
    for (const [start, end] of featureIntervals) {
      for (let i = start; i <= end; ++i) set.add(i);
    }
  }
  return set;
}
