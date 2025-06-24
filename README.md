# PixelFedit

An open-source 8-Bit pixel font editor like I created in 1996.

- Edit a single-byte encoded character set (eg. DOS Codepage 437, 858 or ISO-8859-1)
- supported dimensions are 8x8, 8x16 or 8x16 characters
- save as binary data file, C header or PNG image
- can import png images
- No server-side logic, pure frontend
- Deployed as a GitHub page: <https://learosema.github.io/pixelfedit/>

## Running it local

This project is built using the [Eleventy](https://11ty.dev) Static Site generator.

Run `npm i` && `npm run dev` to install and run Eleventy.

Right now, there is only one single page, but as soon as I feel I want to build more website around this tool, Eleventy will come in handy. Right now, the main use-case of Eleventy is to have a development server with hot refresh. Plus, it comes with easy extensibility to js/css preprocessor toolchains.

