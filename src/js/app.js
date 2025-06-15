import { createSVG, createPNG, createCHeaderFile } from './font-export/index.js';

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

function saveToLocalStorage(state) {
  const data = Array.from(state.font, (v) => v.toString(16).padStart(2, '0')).join('')
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

function loadFont(fileName) {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.responseType = 'arraybuffer';
    xhr.open('GET', fileName);
    xhr.onload = () => {
      if (xhr.status !== 200) {
        reject(xhr.statusText);
        return;
      }
      resolve(new Uint8Array(xhr.response));
    };
    xhr.onerror = (err) => reject(err);
    xhr.send();
  });
}

/**
 * Get a list of fonts that are available on the repo, via GitHub API
 */
function loadOnlineFontList() {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open(
      'GET',
      'https://api.github.com/repos/learosema/pixelfedit/contents/fonts'
    );
    xhr.onload = () => {
      if (xhr.status !== 200) {
        reject(xhr.statusText);
        return;
      }
      try {
        const data = JSON.parse(xhr.response).map((item) => ({
          name: item.name,
          path: item.path,
          cols: 8,
          rows: item.size / 256,
        }));
        resolve(data);
      } catch (ex) {
        reject(ex);
      }
    };
    xhr.onerror = (err) => reject(err);
    xhr.send();
  });
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

function renderCharCursor(state) {
  const pixelSize = state.zoomPixelSize;
  const ctx = charCtx
  ctx.strokeStyle = '#f7c'
  ctx.beginPath()
  ctx.strokeRect(
    state.x * (pixelSize + 1),
    state.y * (pixelSize + 1),
    pixelSize + 2,
    pixelSize + 2
  );
  ctx.stroke()
}

function renderCharacter(state) {
  charCtx.clearRect(0, 0, charCanvas.width, charCanvas.height)

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
    renderCharCursor(state)
  }
}

function togglePixel(state) {
  const line = state.currentChar * state.numRows + state.y;
  state.font[line] ^= 1 << (7 - state.x);
  renderCharacter(state)
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

  loadOnlineFontList().then((data) => {
    const subMenuOpen = $('#subMenuOpen');
    data.map((item) => {
      const $li = $node(`
          <li class="navi__sub-item">
            <a class="navi__sub-link" data-path="${item.path}">${item.name} (8x${item.rows})</a>
          </li>
        `);
      subMenuOpen.appendChild($li);
    });
  });

  const restoredFont = restoreFromLocalStorage()
  if (!restoredFont) {
    loadFont('fonts/8X16.BIN').then(onLoadFont);
  } else {
    onLoadFont(restoredFont);
  }

  $('#menuOpen').addEventListener('click', (e) => {
    $('#menuOpen').classList.toggle('navi__list-item--open');
    if (e.target.hasAttribute('data-path')) {
      loadFont(e.target.getAttribute('data-path')).then(onLoadFont);
    }
    return;
  });

  $('#menuOpenLocal').addEventListener('click', () => {
    const inputFile = $('#inputFile');
    inputFile.click();
  });

  $('#inputFile').addEventListener('change', () => {
    const file = $('#inputFile').files[0];
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

  $('#menuSave').addEventListener('click', () => {
    const anchor = document.createElement('a');
    anchor.setAttribute('download', 'awesome-font.bin');
    anchor.setAttribute(
      'href',
      'data:application/octet-stream;base64,' +
      btoa(String.fromCharCode.apply(null, state.font))
    );
    anchor.click();
  });

  $('#menuExport').addEventListener('click', () => {
    $('#menuExport').classList.toggle('navi__list-item--open');
    return;
  });

  $('#menuExportToSVG').addEventListener('click', (e) => {
    // $('#menuExport').classList.toggle('navi__list-item--open');
    e.preventDefault();
    const anchor = document.createElement('a');
    anchor.setAttribute('download', 'test-svg.svg');
    anchor.setAttribute(
      'href',
      'data:application/octet-stream;base64,' + btoa(createSVG(state.font))
    );
    anchor.click();
  });

  $('#menuExportToTTF').addEventListener('click', (e) => {
    e.preventDefault();
    alert('Temporarily disabled. Needs more work.');
  });

  $('#menuExportToPNG').addEventListener('click', (e) => {
    e.preventDefault();
    const anchor = document.createElement('a');
    anchor.setAttribute('download', 'test-png.png');
    anchor.setAttribute('href', createPNG(state.font));
    anchor.click();
  });

  $('#menuExportToC').addEventListener('click', (e) => {
    e.preventDefault();
    const anchor = document.createElement('a');
    anchor.setAttribute('download', 'fontdata.h');
    anchor.setAttribute('href', createCHeaderFile(state.font));
    anchor.click();
  });


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
    renderFontCanvas(state);
    renderCharacter(state);
  });

  charCanvas.addEventListener('click', (e) => {
    const pixelSize = state.zoomPixelSize;
    const x = e.clientX - charCanvas.offsetLeft + window.scrollX;
    const y = e.clientY - charCanvas.offsetTop + window.scrollY;
    state.x = ((x - 1) / (pixelSize + 1)) | 0;
    state.y = ((y - 1) / (pixelSize + 1)) | 0;
    togglePixel(state);
    renderFontCanvas(state, state.currentChar);
    saveToLocalStorage(state);
  });

  charCanvas.addEventListener('focusin', () => renderCharacter(state))
  charCanvas.addEventListener('focusout', () => renderCharacter(state))

  fontCanvas.addEventListener('keydown', (e) => {
    if (e.code === 'ArrowUp') {
      e.preventDefault();
      state.currentChar-=16;
      if (state.currentChar < 0) {
        state.currentChar += 256
      }
      renderCharacter(state)
      renderFontCanvas(state)
    }
    if (e.code === 'ArrowDown') {
      e.preventDefault();
      state.currentChar += 16;
      if (state.currentChar > 255) {
        state.currentChar -= 256
      }
      renderCharacter(state)
      renderFontCanvas(state)
    }
    if (e.code === 'ArrowLeft') {
      e.preventDefault();
      state.currentChar--;
      if (state.currentChar < 0) {
        state.currentChar += 256
      }
      renderCharacter(state)
      renderFontCanvas(state)
    }
    if (e.code === 'ArrowRight') {
      e.preventDefault();
      state.currentChar++;
      if (state.currentChar > 255) {
        state.currentChar -= 256
      }
      renderCharacter(state)
      renderFontCanvas(state)
    }
  })

  charCanvas.addEventListener('keydown', (e) => {
    if (e.code === 'ArrowUp') {
      e.preventDefault();
      if (state.y > 0) {
        state.y -= 1;
      } else {
        state.y = state.numRows - 1;
      }
      renderCharacter(state)
      return
    }
    if (e.code === 'ArrowDown') {
      e.preventDefault()
      if (state.y < state.numRows - 1) {
        state.y += 1;
      } else {
        state.y = 0;
      }
      renderCharacter(state)
    }

    if (e.code === 'ArrowLeft') {
      e.preventDefault()
      if (state.x > 0) {
        state.x -= 1
      } else {
        state.x = state.numCols - 1
      }
      renderCharacter(state)
      return
    }
    if (e.code === 'ArrowRight') {
      e.preventDefault()
      if (state.x < state.numCols - 1) {
        state.x += 1
      } else {
        state.x = 0
      }
      renderCharacter(state)
    }
    if (e.code === 'Space') {
      togglePixel(state);
      e.preventDefault();
    }
  })
}

setup();
