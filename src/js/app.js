import { createSVG, createPNG, createCHeaderFile } from './font-export/index.js';
import { convertImageData, getImageData, loadImage } from './font-import/read-image.js';
import { generateFromFont } from './generative/gen-from-font.js';

const CP437 = import('./encodings/cp437.js');
const GEN = import('./generative/gen-from-font.js');

const isMac = navigator.userAgent.toLowerCase().includes('mac');
if (isMac) {
  document.body.classList.add('mac');
}
const $ = document.querySelector.bind(document);
const $node = (markup = '<div></div>') => {
  const el = document.createElement('div');
  el.innerHTML = markup.trim();
  const node = el.firstElementChild;
  node.remove();
  return node;
};

const fontCanvas = $('.font-canvas');
const charCanvas = $('.char-canvas');
const fontCtx = fontCanvas.getContext('2d');
const charCtx = charCanvas.getContext('2d');

const state = {
  numRows: 8,
  numCols: 8,
  font: null,
  currentChar: 65,
  pixelSize: 2,
  zoomPixelSize: 24,
  x: 0,
  y: 0,
  onlineFonts: [],
};

function screenreaderAnnounce(message) {
  $('#liveRegion').innerHTML = `<p>${message}</p>`
  console.info('SR: ', message)
}

function screenreaderAnnouncePixel() {
  const bit = 1 << (7 - state.x);
  const line = state.currentChar * state.numRows + state.y;
  screenreaderAnnounce(`${state.x} ${state.y} ${(state.font[line] & bit) > 0 ? 'set':'unset'}`)
}

function clearChar(from, to) {
  if (typeof from !== 'number') {
    from = state.currentChar > -1 ? state.currentChar : 0;
  }
  if (typeof to !== 'number') {
    to = from;
  }
  if (! state.font) {
    return;
  }
  for (let i = from; i <= to; i++) {
    for (let j = 0; j < state.numRows; j++) {
      state.font[i*state.numRows + j] = 0;
    }
  }
  renderFontCanvas();
  renderCharacter();
}

function saveBinary(font, filename = 'awesome-font.bin') {
  const anchor = document.createElement('a');
  anchor.setAttribute('download', filename);
  anchor.setAttribute(
    'href',
    'data:application/octet-stream;base64,' +
    btoa(String.fromCharCode.apply(null, font))
  );
  anchor.click();
}

function savePNG(font, filename = 'awesome-font.png') {
  const anchor = document.createElement('a');
  anchor.setAttribute('download', filename);
  anchor.setAttribute('href', createPNG(font));
  anchor.click();
}

function saveSVG(font, filename = 'awesome-font.svg') {
  const anchor = document.createElement('a');
  anchor.setAttribute('download', filename);
  anchor.setAttribute(
    'href',
    'data:application/octet-stream;base64,' + btoa(createSVG(font))
  );
  anchor.click();
}

function saveCode(font, filename = 'fontdata.h') {
  const anchor = document.createElement('a');
  anchor.setAttribute('download', filename);
  anchor.setAttribute('href', createCHeaderFile(state.font));
  anchor.click();
}

function loadLocal() {
  const inputFile = $('#inputFile');
  inputFile.click();
}

function saveToLocalStorage() {
  const data = Array.from(state.font, (v) => Number(v).toString(16).padStart(2, '0')).join('')
  localStorage.setItem('font', `${state.numCols}x${state.numRows} ${data}`)
}

function restoreFromLocalStorage() {
  const storage = localStorage.getItem('font');
  if (!storage || /^8x\d{1,2} [0-9a-fA-F]{2,16384}$/.test(storage) === false) {
    return null;
  }
  const [, data] = storage.split(' ')
  const font = new Uint8Array(data.length / 2)
  for (let i = 0; i < font.length; i++) {
    font[i] = Number.parseInt(data.slice(i * 2, i * 2 + 2), 16)
  }
  return font
}

const onLoadFont = (font) => {
  state.numRows = (font.length / 256) | 0;
  state.font = font;
  setCanvasSizes();
  renderFontCanvas();
  renderCharacter();
};


/**
 * fetch from url as arraybuffer, reject when http code is 400+
 * @param url
 * @returns {Promise<*>}
 */
async function loadFont(url) {
  if (url.startsWith('data:image') || url.endsWith('.png') || url.endsWith('.svg')) {
    const img = await loadImage(url);
    const imgData = getImageData(img);
    const data = convertImageData(imgData, 8, Math.floor(img.height / 16));
    onLoadFont(data);
    return;
  }

  const response = await fetch(url);
  if (response.status >= 400) {
    throw new Error(`${response.status}: ${response.statusText}`);
  }
  const buf = await response.arrayBuffer();
  onLoadFont(new Uint8ClampedArray(buf));
}

/**
 * Get a list of fonts that are available on the repo, via GitHub API
 */
async function loadOnlineFontList() {
  const response = await fetch(
    'https://api.github.com/repos/learosema/pixelfedit/contents/src/fonts'
  );
  if (response.status >= 400) {
    throw new Error(`${response.status}: ${response.statusText}`);
  }
  return await response.json();
}

/**
 * Render the character table
 * @param {number?} specificChar if specified, only a specific char is updated in the render.
 */
function renderFontCanvas(specificChar) {
  const { pixelSize } = state;
  for (let y = 0; y < 16; y++) {
    for (let x = 0; x < 16; x++) {
      const n = y * 16 + x;
      if (typeof specificChar !== 'undefined' && specificChar !== n) {
        continue;
      }
      fontCtx.fillStyle = state.currentChar === n ? '#f0f' : '#444';
      fontCtx.fillRect(
        x * 10 * pixelSize,
        y * (state.numRows + 2) * pixelSize,
        10 * pixelSize,
        (state.numRows + 2) * pixelSize
      );
      for (let j = 0; j < state.numRows; j++) {
        const line = state.font[n * state.numRows + j];
        for (let i = 0; i < 8; i++) {
          fontCtx.fillStyle = (line & (1 << (7 - i))) > 0 ? '#aaa' : '#000';
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

function renderCharCursor() {
  const pixelSize = state.zoomPixelSize;
  const ctx = charCtx
  ctx.strokeStyle = 'deeppink'
  ctx.lineWidth = 3
  ctx.beginPath()
  ctx.strokeRect(
    state.x * (pixelSize + 1),
    state.y * (pixelSize + 1),
    pixelSize + 2,
    pixelSize + 2
  );
  ctx.stroke()
}

function renderCharacter() {
  charCtx.clearRect(0, 0, charCanvas.width, charCanvas.height)
  CP437.then(({CP437_LABELS}) => {
    $('#characterLabel').textContent =
        `Character ${state.currentChar}: ${CP437_LABELS[state.currentChar]}`;
  })
  const pixelSize = state.zoomPixelSize;
  const n = state.currentChar;
  for (let j = 0; j < state.numRows; j++) {
    const line = state.font[n * state.numRows + j];
    for (let i = 0; i < 8; i++) {
      charCtx.fillStyle = (line & (1 << (7 - i))) > 0 ? '#aaa' : '#000';
      charCtx.fillRect(
        1 + i * (pixelSize + 1),
        1 + j * (pixelSize + 1),
        pixelSize,
        pixelSize
      );
    }
  }
  if (document.activeElement === charCanvas) {
    renderCharCursor()
  }
}

function togglePixel() {
  const bit = 1 << (7 - state.x);
  const line = state.currentChar * state.numRows + state.y;
  state.font[line] ^= bit;
  renderCharacter()
  screenreaderAnnouncePixel()
}

function setCanvasSizes() {
  fontCanvas.width = 10 * state.pixelSize * 16;
  fontCanvas.height = (state.numRows + 2) * state.pixelSize * 16;
  charCanvas.width = 1 + 8 * (state.zoomPixelSize + 1);
  charCanvas.height = 1 + state.numRows * (state.zoomPixelSize + 1);
}

function populateFontList(entries) {
  const selectOpen = $('#selectOpen');
  selectOpen.innerHTML = '<option value="local">local file</option>';
  for (const item of entries) {
    const $li = $node(`
      <option value="fonts/${item.name}">${item.name}</option>
    `);
    selectOpen.appendChild($li);
  }
}

function initCharCanvasEvents() {
  charCanvas.addEventListener('click', (e) => {
    const pixelSize = state.zoomPixelSize;
    const x = e.clientX - charCanvas.offsetLeft + window.scrollX;
    const y = e.clientY - charCanvas.offsetTop + window.scrollY;
    state.x = ((x - 1) / (pixelSize + 1)) | 0;
    state.y = ((y - 1) / (pixelSize + 1)) | 0;
    togglePixel();
    renderCharacter();
    renderFontCanvas(state.currentChar);
    saveToLocalStorage();
  });

  charCanvas.addEventListener('focusin', () => renderCharacter());
  charCanvas.addEventListener('focusout', () => renderCharacter());

  charCanvas.addEventListener('keydown', (e) => {
    if (e.code === 'ArrowUp') {
      e.preventDefault();
      if (state.y > 0) {
        state.y -= 1;
      } else {
        state.y = state.numRows - 1;
      }
      renderCharacter()
      screenreaderAnnouncePixel()
      return
    }
    if (e.code === 'ArrowDown') {
      e.preventDefault()
      if (state.y < state.numRows - 1) {
        state.y += 1;
      } else {
        state.y = 0;
      }
      renderCharacter()
      screenreaderAnnouncePixel()
    }

    if (e.code === 'ArrowLeft') {
      e.preventDefault()
      if (state.x > 0) {
        state.x -= 1
      } else {
        state.x = state.numCols - 1
      }
      renderCharacter()
      screenreaderAnnouncePixel()
    }
    if (e.code === 'ArrowRight') {
      e.preventDefault()
      if (state.x < state.numCols - 1) {
        state.x += 1
      } else {
        state.x = 0
      }
      renderCharacter()
      screenreaderAnnouncePixel()
    }
    if (e.code === 'Space') {
      togglePixel();
      renderFontCanvas();
      e.preventDefault();
    }
  })
}

function openCommandPalette() {
  /** @type HTMLInputElement */
  const input = $('#commandInput');
  $('#commandPaletteDialog').showModal();
  input.value = '';
  input.focus();
}

function initMenu() {
  $('#menuOpen').addEventListener('click', () => {
    loadOnlineFontList().then(populateFontList);
    $('#openDialog').showModal();
  });

  $('#menuCommand').addEventListener('click', () => {
    openCommandPalette();
  })

  $('#menuSave').addEventListener('click', () => {
    $('#saveDialog').showModal();
  });
  import('tinykeys').then(({tinykeys}) => {
    tinykeys(window, {"$mod+P": event => {
      event.preventDefault();
        openCommandPalette();
      }})
  })
  $('#commandForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const cmd = $('#commandInput').value;
    if (cmd === 'open') {
      setTimeout(() => $('#openDialog').showModal(), 0)
    }
    if (cmd === 'save') {
      setTimeout(() => $('#saveDialog').showModal(), 0)
    }
    if (cmd.startsWith('open ')) {
      const fontUrl = cmd.slice(5);
      if (fontUrl.startsWith('https://') || fontUrl.startsWith('http://')) {
        loadFont(fontUrl).then(() => console.info('Font loaded'));
      } else {
        loadFont(`fonts/${fontUrl}`).then(() =>
          console.info('Font loaded'));
      }
    }
    if (cmd === 'clear') {
       clearChar(state.currentChar);
    }
    if (cmd.startsWith('clear ')) {
      const param = cmd.slice(6).split(/[- ]/)
        .map(str => str.trim())
        .filter(str => str !== '');
      if (param.length === 3 && param[0] === 'except') {
        clearChar(0, (Number(param[1]) || 1) - 1);
        clearChar((Number(param[2]) || 0) + 1, 255);
      }
      if (param.length === 1) {
        clearChar(Number(param[0]) || 0);
      }
      if (param.length === 2) {
        clearChar(Number(param[0]) || 0, Number(param[1]) || 0);
      }
    }
    if (cmd.startsWith('serif ')) {
      const param = cmd.slice(6).split(/[- ]/)
        .map(str => str.trim())
        .filter(str => str !== '');
      if (param.length === 3 && param[0] === 'except') {
        const excludeA =  (Number(param[1]) || 1)
        const excludeB =  (Number(param[2]) || 1)

        generateFromFont(state,
          'serif',
          0, excludeA - 1,
          8, state.numRows);

        generateFromFont(state,
          'serif',
          excludeB + 1, 255,
          8, state.numRows);
      }
      if (param.length === 1) {
        const _to = Number(param[0]) || 0
        generateFromFont(state, 'serif', _to, _to);
      }
      if (param.length === 2) {
        const from = Number(param[0]) || 0
        const _to = Number(param[1]) || 0
        generateFromFont(state, 'serif', from, _to);
      }
      renderCharacter();
      renderFontCanvas();
    }
    if (cmd.startsWith('sans-serif ')) {
      const param = cmd.slice(11).split(/[- ]/)
        .map(str => str.trim())
        .filter(str => str !== '');
      if (param.length === 3 && param[0] === 'except') {
        const excludeA =  (Number(param[1]) || 1)
        const excludeB =  (Number(param[2]) || 1)

        generateFromFont(state,
          'sans-serif',
          0, excludeA - 1,
          8, state.numRows);

        generateFromFont(state,
          'sans-serif',
          excludeB + 1, 255,
          8, state.numRows);
      }
      if (param.length === 1) {
        const _to = Number(param[0]) || 0
        generateFromFont(state, 'sans-serif', _to, _to);
      }
      if (param.length === 2) {
        const from = Number(param[0]) || 0
        const _to = Number(param[1]) || 0
        generateFromFont(state, 'sans-serif', from, _to);
      }
      renderCharacter();
      renderFontCanvas();
    }




    $('#commandPaletteDialog').close();
  })
}

function initFontCanvasEvents() {
  fontCanvas.addEventListener('click', (e) => {
    if (!state.font) {
      return;
    }
    const { pixelSize, numRows } = state;
    const x = e.clientX - fontCanvas.offsetLeft + window.scrollX;
    const y = e.clientY - fontCanvas.offsetTop + window.scrollY;
    const xIndex = (x / (10 * pixelSize)) | 0;
    const yIndex = (y / ((numRows + 2) * pixelSize)) | 0;
    state.currentChar = yIndex * 16 + xIndex;
    renderFontCanvas();
    renderCharacter();
  });

  fontCanvas.addEventListener('keydown', (e) => {
    if (e.code === 'ArrowUp') {
      e.preventDefault();
      state.currentChar -= 16;
      if (state.currentChar < 0) {
        state.currentChar += 256;
      }
      renderCharacter();
      renderFontCanvas();
    }
    if (e.code === 'ArrowDown') {
      e.preventDefault();
      state.currentChar += 16;
      if (state.currentChar > 255) {
        state.currentChar -= 256;
      }
      renderCharacter();
      renderFontCanvas();
    }
    if (e.code === 'ArrowLeft') {
      e.preventDefault();
      state.currentChar--;
      if (state.currentChar < 0) {
        state.currentChar += 256;
      }
      renderCharacter();
      renderFontCanvas();
    }
    if (e.code === 'ArrowRight') {
      e.preventDefault();
      state.currentChar++;
      if (state.currentChar > 255) {
        state.currentChar -= 256;
      }
      renderCharacter();
      renderFontCanvas();
    }
  });
}

function initOpenDialog() {
  $('#openForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const selectOpen = $('#selectOpen').value;
    if (selectOpen === 'local') {
      loadLocal();
    } else {
      loadFont(selectOpen);
    }
    $('#openDialog').close();
  });

  $('#inputFile').addEventListener('change', () => {
    const file = $('#inputFile').files[0];
    if (!file || !(file instanceof Blob)) {
      return;
    }
    if (file.name.endsWith('.png')) {
      const reader = new FileReader();
      reader.onload = async () => {
        await loadFont(reader.result);
      };
      reader.readAsDataURL(file);
      return;
    }
    if (file.name.endsWith('.bin')) {
      const reader = new FileReader();
      reader.onload = () => {
        onLoadFont(new Uint8Array(reader.result));
      };
      reader.readAsArrayBuffer(file);
    }
  });
}

function initSaveDialog() {
  $('#saveForm').addEventListener('submit', (e) => {
    e.preventDefault();
    const selectSaveFormat = $('#selectSaveFormat').value;
    switch (selectSaveFormat) {
      case 'binary':
        saveBinary(state.font);
        break;
      case 'png':
        savePNG(state.font);
        break;
      case 'svg':
        saveSVG(state.font);
        break;
      case 'c':
        saveCode(state.font);
    }
    $('#saveDialog').close();
  });
}

async function setup() {
  const restoredFont = restoreFromLocalStorage()
  if (!restoredFont) {
    await loadFont('fonts/8X16.BIN');
  } else {
    onLoadFont(restoredFont);
  }
  initCharCanvasEvents();
  initMenu();
  initSaveDialog();
  initOpenDialog();
  initFontCanvasEvents();
}

setup().then(() => console.info('initialization successful'));
