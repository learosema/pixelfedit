export class GlyphEditor extends HTMLElement {

  svg = null;

  static register() {
    customElements.define('glyph-editor', GlyphEditor);
  }

  constructor() {
    super();
    this.style.display = 'block';
  }

  connectedCallback() {
    
    if (! this.svg) {
      this.svg = this.querySelector('svg.glyph__canvas')
    }
    if (! this.svg) {
      this.svg = document.createElement('svg')
      this.svg.classList.add('glyph__canvas')
      this.appendChild(this.svg)
    }
  }

  

}
