const $ = document.querySelector.bind(document);
const fontCanvas = $(".font-canvas");
const charCanvas = $(".char-canvas");
const fontCtx = fontCanvas.getContext("2d");
const charCtx = charCanvas.getContext("2d");

const state = {
  numRows: 8,
  font: null,
  currentChar: 65,
  pixelSize: 2,
  zoomPixelSize: 24,
  x: 0,
  y: 0,
};

function loadFont(fileName) {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.responseType = "arraybuffer";
    xhr.open("GET", fileName);
    xhr.onload = () => {
      resolve(new Uint8Array(xhr.response));
    };
    xhr.onerror = (err) => reject(err);
    xhr.send();
  });
}

function hexDump(data) {
  return Array.prototype.concat
    .apply([], data)
    .map((c) => c.toString(16))
    .join("");
}

/**
 * Render the character table
 * @param {*} state the state object
 * @param {*} specificChar if specified, only a specific char is updated in the render.
 */
function renderFontCanvas(state, specificChar) {
  const { pixelSize } = state;
  for (let y = 0; y < 16; y++) {
    for (let x = 0; x < 16; x++) {
      const n = y * 16 + x;
      if (typeof specificChar !== "undefined" && specificChar !== n) {
        continue;
      }
      fontCtx.fillStyle = state.currentChar === n ? "#f0f" : "#444";
      fontCtx.fillRect(
        x * 10 * pixelSize,
        y * (state.numRows + 2) * pixelSize,
        10 * pixelSize,
        (state.numRows + 2) * pixelSize
      );
      for (let j = 0; j < state.numRows; j++) {
        const line = state.font[n * state.numRows + j];
        for (let i = 0; i < 8; i++) {
          fontCtx.fillStyle = (line & (1 << (7 - i))) > 0 ? "#aaa" : "#000";
          fontCtx.fillRect(
            (x * 10 + i + 1) * pixelSize,
            (y * (state.numRows + 2) + j + 1) * pixelSize,
            pixelSize,
            pixelSize
          );
        }
      }
    }
  }
}

function renderCharacter(state) {
  const pixelSize = state.zoomPixelSize;
  const n = state.currentChar;
  for (let j = 0; j < state.numRows; j++) {
    const line = state.font[n * state.numRows + j];
    for (let i = 0; i < 8; i++) {
      charCtx.fillStyle = (line & (1 << (7 - i))) > 0 ? "#aaa" : "#000";
      charCtx.fillRect(
        1 + i * (pixelSize + 1),
        1 + j * (pixelSize + 1),
        pixelSize,
        pixelSize
      );
    }
  }
}

function setCanvasSizes(state) {
  fontCanvas.width = 10 * state.pixelSize * 16;
  fontCanvas.height = (state.numRows + 2) * state.pixelSize * 16;
  charCanvas.width = 1 + 8 * (state.zoomPixelSize + 1);
  charCanvas.height = 1 + state.numRows * (state.zoomPixelSize + 1);
}

function setup() {
  const onLoadFont = (font) => {
    state.numRows = (font.length / 256) | 0;
    state.font = font;
    setCanvasSizes(state);
    renderFontCanvas(state);
    renderCharacter(state);
  };

  loadFont("8X16.BIN").then(onLoadFont);

  $("#menuOpen").addEventListener("click", (e) => {
    console.log("muh", e, e.target);
    const inputFile = $("#inputFile");
    inputFile.click();
  });

  $("#inputFile").addEventListener("change", (e) => {
    const file = $("#inputFile").files[0];
    // TODO: if file.name ends with .png then try to import?
    if (!file || file instanceof Blob === false) {
      return;
    }
    const reader = new FileReader();
    reader.onload = () => {
      console.log(reader.result);
      onLoadFont(new Uint8Array(reader.result));
    };
    reader.readAsArrayBuffer(file);
  });

  fontCanvas.addEventListener("click", (e) => {
    if (!state.font) {
      return;
    }
    const { pixelSize, numRows } = state;
    const x = e.clientX - fontCanvas.offsetLeft + window.scrollX;
    const y = e.clientY - fontCanvas.offsetTop + window.scrollY;
    const xIndex = (x / (10 * pixelSize)) | 0;
    const yIndex = (y / ((numRows + 2) * pixelSize)) | 0;
    state.currentChar = yIndex * 16 + xIndex;
    renderFontCanvas(state);
    renderCharacter(state);
  });

  charCanvas.addEventListener("click", (e) => {
    const pixelSize = state.zoomPixelSize;
    const x = e.clientX - charCanvas.offsetLeft + window.scrollX;
    const y = e.clientY - charCanvas.offsetTop + window.scrollY;
    const xIndex = ((x - 1) / (pixelSize + 1)) | 0;
    const yIndex = ((y - 1) / (pixelSize + 1)) | 0;
    console.log(xIndex, yIndex);
    const line = state.currentChar * state.numRows + yIndex;
    state.font[line] ^= 1 << (7 - xIndex);
    renderCharacter(state);
    renderFontCanvas(state, state.currentChar);
  });
}

setup();
